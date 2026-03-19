#!/usr/bin/env python3

import argparse
import csv
import json
import math
import os
import random
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path


CHILD_CODE = r"""
import json
import sys
import time
from pathlib import Path

module_dir = sys.argv[1]
task = sys.argv[2]
payload_path = Path(sys.argv[3])
sys.path.insert(0, module_dir)
import RNA

payload = json.loads(payload_path.read_text())
workload = payload["workload"]
repeats = payload["repeats"]
warmups = payload["warmups"]
results = []

for item in workload:
    identifier = item["id"]
    durations = []
    outputs = None
    for _ in range(warmups):
        if task == "fold":
            RNA.fold(item["sequence"])
        elif task == "cofold":
            RNA.cofold(item["sequence"])
        elif task == "subopt":
            RNA.subopt(item["sequence"], item["delta"])
    for _ in range(repeats):
        t0 = time.perf_counter_ns()
        if task == "fold":
            structure, energy = RNA.fold(item["sequence"])
            current = {
                "structure": structure,
                "energy": float(energy),
            }
        elif task == "cofold":
            structure, energy = RNA.cofold(item["sequence"])
            current = {
                "structure": structure,
                "energy": float(energy),
            }
        elif task == "subopt":
            sols = RNA.subopt(item["sequence"], item["delta"])
            current = {
                "solutions": [
                    {
                        "structure": s.structure,
                        "energy": float(s.energy),
                    } for s in sols
                ]
            }
        else:
            raise RuntimeError(f"unknown task: {task}")
        durations.append(time.perf_counter_ns() - t0)
        outputs = current
    results.append({
        "id": identifier,
        "meta": item,
        "durations_ns": durations,
        "output": outputs,
    })

print(json.dumps(results))
"""


def synthetic_sequence(length: int, seed: int) -> str:
    alphabet = "ACGU"
    state = seed ^ length
    out = []
    for _ in range(length):
        state = (state * 1664525 + 1013904223) & 0xFFFFFFFF
        out.append(alphabet[(state >> 24) & 0x3])
    return "".join(out)


def fold_workload():
    workload = []
    lengths = [100, 500, 1000]
    idx = 1
    for length in lengths:
        for rep in range(30):
            workload.append({
                "id": f"fold_{idx:03d}",
                "length": length,
                "sequence": synthetic_sequence(length, 1729 + rep * 17 + length),
            })
            idx += 1
    return workload


def cofold_workload():
    workload = []
    lengths = [50, 100, 200]
    pairs = [(a, b) for a in lengths for b in lengths]
    idx = 1
    for rep in range(30):
        l1, l2 = pairs[rep % len(pairs)]
        s1 = synthetic_sequence(l1, 1729 + rep * 31 + l1)
        s2 = synthetic_sequence(l2, 3141 + rep * 47 + l2)
        workload.append({
            "id": f"cofold_{idx:03d}",
            "len1": l1,
            "len2": l2,
            "sequence": f"{s1}&{s2}",
        })
        idx += 1
    return workload


def subopt_workload():
    workload = []
    lengths = [50, 100, 150, 200]
    deltas = [5, 10]
    idx = 1
    for rep in range(30):
        length = lengths[rep % len(lengths)]
        delta = deltas[(rep // len(lengths)) % len(deltas)]
        workload.append({
            "id": f"subopt_{idx:03d}",
            "length": length,
            "delta": delta,
            "sequence": synthetic_sequence(length, 2718 + rep * 29 + length),
        })
        idx += 1
    return workload


def run_module(module_dir: Path, task: str, workload, repeats: int, warmups: int):
    payload = {
        "workload": workload,
        "repeats": repeats,
        "warmups": warmups,
    }
    with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as tf:
        tf.write(json.dumps(payload))
        payload_path = Path(tf.name)
    try:
        proc = subprocess.run(
            [sys.executable, "-c", CHILD_CODE, str(module_dir), task, str(payload_path)],
            check=True,
            capture_output=True,
            text=True,
        )
        return json.loads(proc.stdout)
    finally:
        payload_path.unlink(missing_ok=True)


def compare_fold_or_cofold(old_rows, new_rows):
    diffs = []
    csv_rows = []
    for old, new in zip(old_rows, new_rows):
        same_struct = old["output"]["structure"] == new["output"]["structure"]
        energy_delta = new["output"]["energy"] - old["output"]["energy"]
        same_energy = abs(energy_delta) < 1e-9
        if not (same_struct and same_energy):
            diffs.append(old["id"])
        csv_rows.append({
            "id": old["id"],
            **{k: v for k, v in old["meta"].items() if k not in {"id", "sequence"}},
            "sequence": old["meta"]["sequence"],
            "old_structure": old["output"]["structure"],
            "old_energy": old["output"]["energy"],
            "new_structure": new["output"]["structure"],
            "new_energy": new["output"]["energy"],
            "structure_match": int(same_struct),
            "energy_delta": energy_delta,
            "old_median_ns": int(statistics.median(old["durations_ns"])),
            "new_median_ns": int(statistics.median(new["durations_ns"])),
        })
    return diffs, csv_rows


def compare_subopt(old_rows, new_rows):
    diffs = []
    csv_rows = []
    for old, new in zip(old_rows, new_rows):
        old_solutions = old["output"]["solutions"]
        new_solutions = new["output"]["solutions"]
        same = old_solutions == new_solutions
        if not same:
            diffs.append(old["id"])
        csv_rows.append({
            "id": old["id"],
            "length": old["meta"]["length"],
            "delta": old["meta"]["delta"],
            "sequence": old["meta"]["sequence"],
            "old_solution_count": len(old_solutions),
            "new_solution_count": len(new_solutions),
            "solutions_match": int(same),
            "old_first_structure": old_solutions[0]["structure"] if old_solutions else "",
            "old_first_energy": old_solutions[0]["energy"] if old_solutions else "",
            "new_first_structure": new_solutions[0]["structure"] if new_solutions else "",
            "new_first_energy": new_solutions[0]["energy"] if new_solutions else "",
            "old_median_ns": int(statistics.median(old["durations_ns"])),
            "new_median_ns": int(statistics.median(new["durations_ns"])),
        })
    return diffs, csv_rows


def write_csv(path: Path, rows):
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def summarize(task, old_rows, new_rows, diffs, bootstrap_iterations, significance_threshold):
    old_medians = [int(statistics.median(r["durations_ns"])) for r in old_rows]
    new_medians = [int(statistics.median(r["durations_ns"])) for r in new_rows]
    old_total = sum(old_medians)
    new_total = sum(new_medians)
    ratio_samples = bootstrap_runtime_ratio(old_medians, new_medians, iterations=bootstrap_iterations)
    ci_low, ci_high = percentile_bounds(ratio_samples, 0.95)
    ratio = (new_total / old_total) if old_total else None
    return {
        "task": task,
        "cases": len(old_rows),
        "differences": len(diffs),
        "old_total_median_ns": old_total,
        "new_total_median_ns": new_total,
        "runtime_ratio": ratio,
        "runtime_ratio_ci95_low": ci_low,
        "runtime_ratio_ci95_high": ci_high,
        "significance_threshold": significance_threshold,
        "speedup_significant": bool(ci_high < significance_threshold) if ci_high is not None else False,
        "difference_ids": diffs,
    }


def bootstrap_runtime_ratio(old_medians, new_medians, iterations=5000, seed=1729):
    if not old_medians or not new_medians or len(old_medians) != len(new_medians):
      return []
    rng = random.Random(seed)
    ratios = []
    count = len(old_medians)
    for _ in range(iterations):
        sample_old = 0
        sample_new = 0
        for _ in range(count):
            idx = rng.randrange(count)
            sample_old += old_medians[idx]
            sample_new += new_medians[idx]
        ratios.append(sample_new / sample_old if sample_old else math.nan)
    return [r for r in ratios if math.isfinite(r)]


def percentile_bounds(values, confidence):
    if not values:
        return None, None
    sorted_values = sorted(values)
    low_q = (1.0 - confidence) / 2.0
    high_q = 1.0 - low_q
    low_idx = max(0, min(len(sorted_values) - 1, int(low_q * (len(sorted_values) - 1))))
    high_idx = max(0, min(len(sorted_values) - 1, int(high_q * (len(sorted_values) - 1))))
    return sorted_values[low_idx], sorted_values[high_idx]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--new-module-dir", required=True)
    parser.add_argument("--old-module-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmups", type=int, default=0)
    parser.add_argument("--bootstrap-iterations", type=int, default=5000)
    parser.add_argument("--significance-threshold", type=float, default=0.95)
    parser.add_argument("--tasks", default="fold,cofold,subopt")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    all_tasks = {
        "fold": fold_workload(),
        "cofold": cofold_workload(),
        "subopt": subopt_workload(),
    }
    requested_tasks = [task.strip() for task in args.tasks.split(",") if task.strip()]
    tasks = {task: all_tasks[task] for task in requested_tasks}

    overall = {}
    for task, workload in tasks.items():
        old_rows = run_module(Path(args.old_module_dir), task, workload, args.repeats, args.warmups)
        new_rows = run_module(Path(args.new_module_dir), task, workload, args.repeats, args.warmups)

        if task in {"fold", "cofold"}:
            diffs, csv_rows = compare_fold_or_cofold(old_rows, new_rows)
        else:
            diffs, csv_rows = compare_subopt(old_rows, new_rows)

        write_csv(out_dir / f"python_{task}_compare.csv", csv_rows)
        summary = summarize(task,
                            old_rows,
                            new_rows,
                            diffs,
                            args.bootstrap_iterations,
                            args.significance_threshold)
        overall[task] = summary
        (out_dir / f"python_{task}_compare.summary.json").write_text(json.dumps(summary, indent=2))

    (out_dir / "python_viennarna_compare.summary.json").write_text(json.dumps(overall, indent=2))
    print(json.dumps(overall, indent=2))


if __name__ == "__main__":
    main()

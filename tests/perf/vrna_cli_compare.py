#!/usr/bin/env python3

import argparse
import csv
import json
import math
import random
import re
import statistics
import subprocess
import time
from pathlib import Path


RESULT_RE = re.compile(r"^([().&]+)\s+\(\s*([-+]?\d+(?:\.\d+)?)\)")
SUBOPT_RE = re.compile(r"^([().&]+)\s+([-+]?\d+(?:\.\d+)?)$")


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
    idx = 1
    for length in (100, 500, 1000):
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
    lengths = (50, 100, 200)
    pairs = [(a, b) for a in lengths for b in lengths]
    for idx in range(30):
        l1, l2 = pairs[idx % len(pairs)]
        workload.append({
            "id": f"cofold_{idx + 1:03d}",
            "len1": l1,
            "len2": l2,
            "sequence1": synthetic_sequence(l1, 1729 + idx * 31 + l1),
            "sequence2": synthetic_sequence(l2, 3141 + idx * 47 + l2),
        })
    return workload


def subopt_workload():
    workload = []
    lengths = (50, 100, 150, 200)
    deltas = (5, 10)
    for idx in range(30):
        length = lengths[idx % len(lengths)]
        delta = deltas[(idx // len(lengths)) % len(deltas)]
        workload.append({
            "id": f"subopt_{idx + 1:03d}",
            "length": length,
            "delta": delta,
            "sequence": synthetic_sequence(length, 2718 + idx * 29 + length),
        })
    return workload


def run_fold(binary: Path, sequence: str):
    start = time.perf_counter_ns()
    proc = subprocess.run([str(binary), "--noPS"],
                          input=f"{sequence}\n",
                          text=True,
                          capture_output=True,
                          check=True)
    elapsed = time.perf_counter_ns() - start
    lines = [line.strip() for line in proc.stdout.splitlines() if line.strip()]
    match = RESULT_RE.match(lines[-1])
    return {
        "structure": match.group(1),
        "energy": float(match.group(2)),
    }, elapsed


def run_cofold(binary: Path, sequence1: str, sequence2: str):
    start = time.perf_counter_ns()
    proc = subprocess.run([str(binary), "--noPS"],
                          input=f"{sequence1}&{sequence2}\n",
                          text=True,
                          capture_output=True,
                          check=True)
    elapsed = time.perf_counter_ns() - start
    lines = [line.strip() for line in proc.stdout.splitlines() if line.strip()]
    match = RESULT_RE.match(lines[-1])
    return {
        "structure": match.group(1),
        "energy": float(match.group(2)),
    }, elapsed


def run_subopt(binary: Path, sequence: str, delta: int):
    start = time.perf_counter_ns()
    proc = subprocess.run([str(binary), "-e", f"{delta / 100.0:.2f}"],
                          input=f"{sequence}\n",
                          text=True,
                          capture_output=True,
                          check=True)
    elapsed = time.perf_counter_ns() - start
    lines = [line.strip() for line in proc.stdout.splitlines() if line.strip()]
    solutions = []
    for line in lines[1:]:
        match = SUBOPT_RE.match(line)
        if match:
            solutions.append({
                "structure": match.group(1),
                "energy": float(match.group(2)),
            })
    return {"solutions": solutions}, elapsed


def run_task(binary: Path, task: str, workload, repeats: int, warmups: int):
    rows = []
    for item in workload:
        for _ in range(warmups):
            if task == "fold":
                run_fold(binary, item["sequence"])
            elif task == "cofold":
                run_cofold(binary, item["sequence1"], item["sequence2"])
            else:
                run_subopt(binary, item["sequence"], item["delta"])

        durations = []
        output = None
        for _ in range(repeats):
            if task == "fold":
                output, elapsed = run_fold(binary, item["sequence"])
                meta = {
                    "length": item["length"],
                    "sequence": item["sequence"],
                }
            elif task == "cofold":
                output, elapsed = run_cofold(binary, item["sequence1"], item["sequence2"])
                meta = {
                    "len1": item["len1"],
                    "len2": item["len2"],
                    "sequence1": item["sequence1"],
                    "sequence2": item["sequence2"],
                }
            else:
                output, elapsed = run_subopt(binary, item["sequence"], item["delta"])
                meta = {
                    "length": item["length"],
                    "delta": item["delta"],
                    "sequence": item["sequence"],
                }
            durations.append(elapsed)
        rows.append({
            "id": item["id"],
            "meta": meta,
            "durations_ns": durations,
            "output": output,
        })
    return rows


def bootstrap_runtime_ratio(old_medians, new_medians, iterations=5000, seed=1729):
    if not old_medians or len(old_medians) != len(new_medians):
        return []
    rng = random.Random(seed)
    ratios = []
    count = len(old_medians)
    for _ in range(iterations):
        so = 0
        sn = 0
        for _ in range(count):
            idx = rng.randrange(count)
            so += old_medians[idx]
            sn += new_medians[idx]
        ratios.append(sn / so if so else math.nan)
    return [r for r in ratios if math.isfinite(r)]


def percentile_bounds(values, confidence):
    if not values:
        return None, None
    values = sorted(values)
    low_idx = int(((1.0 - confidence) / 2.0) * (len(values) - 1))
    high_idx = int((1.0 - ((1.0 - confidence) / 2.0)) * (len(values) - 1))
    return values[low_idx], values[high_idx]


def compare_fold_like(old_rows, new_rows):
    diffs = []
    csv_rows = []
    for old, new in zip(old_rows, new_rows):
        same_struct = old["output"]["structure"] == new["output"]["structure"]
        delta = new["output"]["energy"] - old["output"]["energy"]
        same_energy = abs(delta) < 1e-9
        if not (same_struct and same_energy):
            diffs.append(old["id"])
        csv_rows.append({
            "id": old["id"],
            **old["meta"],
            "old_structure": old["output"]["structure"],
            "old_energy": old["output"]["energy"],
            "new_structure": new["output"]["structure"],
            "new_energy": new["output"]["energy"],
            "structure_match": int(same_struct),
            "energy_delta": delta,
            "old_median_ns": int(statistics.median(old["durations_ns"])),
            "new_median_ns": int(statistics.median(new["durations_ns"])),
        })
    return diffs, csv_rows


def compare_subopt(old_rows, new_rows):
    diffs = []
    csv_rows = []
    for old, new in zip(old_rows, new_rows):
        same = old["output"]["solutions"] == new["output"]["solutions"]
        if not same:
            diffs.append(old["id"])
        csv_rows.append({
            "id": old["id"],
            **old["meta"],
            "old_solution_count": len(old["output"]["solutions"]),
            "new_solution_count": len(new["output"]["solutions"]),
            "solutions_match": int(same),
            "old_first_structure": old["output"]["solutions"][0]["structure"] if old["output"]["solutions"] else "",
            "old_first_energy": old["output"]["solutions"][0]["energy"] if old["output"]["solutions"] else "",
            "new_first_structure": new["output"]["solutions"][0]["structure"] if new["output"]["solutions"] else "",
            "new_first_energy": new["output"]["solutions"][0]["energy"] if new["output"]["solutions"] else "",
            "old_median_ns": int(statistics.median(old["durations_ns"])),
            "new_median_ns": int(statistics.median(new["durations_ns"])),
        })
    return diffs, csv_rows


def summarize(task, old_rows, new_rows, diffs, bootstrap_iterations, significance_threshold):
    old_medians = [int(statistics.median(r["durations_ns"])) for r in old_rows]
    new_medians = [int(statistics.median(r["durations_ns"])) for r in new_rows]
    ratio_samples = bootstrap_runtime_ratio(old_medians, new_medians, iterations=bootstrap_iterations)
    ci_low, ci_high = percentile_bounds(ratio_samples, 0.95)
    old_total = sum(old_medians)
    new_total = sum(new_medians)
    return {
        "task": task,
        "cases": len(old_rows),
        "differences": len(diffs),
        "old_total_median_ns": old_total,
        "new_total_median_ns": new_total,
        "runtime_ratio": (new_total / old_total) if old_total else None,
        "runtime_ratio_ci95_low": ci_low,
        "runtime_ratio_ci95_high": ci_high,
        "significance_threshold": significance_threshold,
        "speedup_significant": bool(ci_high < significance_threshold) if ci_high is not None else False,
        "difference_ids": diffs,
    }


def write_csv(path: Path, rows):
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--new-rnafold-bin", required=True)
    parser.add_argument("--old-rnafold-bin", required=True)
    parser.add_argument("--new-rnacofold-bin", required=True)
    parser.add_argument("--old-rnacofold-bin", required=True)
    parser.add_argument("--new-rnasubopt-bin", required=True)
    parser.add_argument("--old-rnasubopt-bin", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--tasks", default="fold,cofold,subopt")
    parser.add_argument("--warmups", type=int, default=0)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--bootstrap-iterations", type=int, default=5000)
    parser.add_argument("--significance-threshold", type=float, default=0.95)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    workloads = {
        "fold": fold_workload(),
        "cofold": cofold_workload(),
        "subopt": subopt_workload(),
    }
    binaries = {
        "fold": (Path(args.old_rnafold_bin), Path(args.new_rnafold_bin)),
        "cofold": (Path(args.old_rnacofold_bin), Path(args.new_rnacofold_bin)),
        "subopt": (Path(args.old_rnasubopt_bin), Path(args.new_rnasubopt_bin)),
    }
    tasks = [task.strip() for task in args.tasks.split(",") if task.strip()]

    overall = {}
    for task in tasks:
        old_rows = run_task(binaries[task][0], task, workloads[task], args.repeats, args.warmups)
        new_rows = run_task(binaries[task][1], task, workloads[task], args.repeats, args.warmups)
        if task == "subopt":
            diffs, rows = compare_subopt(old_rows, new_rows)
        else:
            diffs, rows = compare_fold_like(old_rows, new_rows)
        write_csv(out_dir / f"cli_{task}_compare.csv", rows)
        summary = summarize(task,
                            old_rows,
                            new_rows,
                            diffs,
                            args.bootstrap_iterations,
                            args.significance_threshold)
        (out_dir / f"cli_{task}_compare.summary.json").write_text(json.dumps(summary, indent=2),
                                                                  encoding="utf-8")
        overall[task] = summary

    (out_dir / "cli_compare.summary.json").write_text(json.dumps(overall, indent=2), encoding="utf-8")
    print(json.dumps(overall, indent=2))


if __name__ == "__main__":
    main()

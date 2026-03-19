#!/usr/bin/env python3

import argparse
import csv
import json
import math
import random
import statistics
import subprocess
import tempfile
from pathlib import Path


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
        s1 = synthetic_sequence(l1, 1729 + idx * 31 + l1)
        s2 = synthetic_sequence(l2, 3141 + idx * 47 + l2)
        workload.append({
            "id": f"cofold_{idx + 1:03d}",
            "len1": l1,
            "len2": l2,
            "sequence": f"{s1}&{s2}",
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


def write_workload(path: Path, task: str, workload):
    lines = []
    for item in workload:
        if task == "subopt":
            lines.append(f"{item['id']}\t{item['delta']}\t{item['sequence']}")
        else:
            lines.append(f"{item['id']}\t{item['sequence']}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_worker(worker: Path, task: str, workload, repeats: int, warmups: int):
    with tempfile.TemporaryDirectory(prefix="vrna-capi-") as tmpdir:
        input_path = Path(tmpdir) / "input.tsv"
        output_path = Path(tmpdir) / "output.json"
        write_workload(input_path, task, workload)
        subprocess.run(
            [str(worker), task, str(input_path), str(output_path), str(repeats), str(warmups)],
            check=True,
        )
        return json.loads(output_path.read_text(encoding="utf-8"))["results"]


def run_probe(worker: Path):
    with tempfile.TemporaryDirectory(prefix="vrna-probe-") as tmpdir:
        output_path = Path(tmpdir) / "probe.json"
        subprocess.run([str(worker), "probe", str(output_path)], check=True)
        return json.loads(output_path.read_text(encoding="utf-8"))


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
          **{k: v for k, v in old["meta"].items()},
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
            **{k: v for k, v in old["meta"].items()},
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
    parser.add_argument("--new-worker", required=True)
    parser.add_argument("--old-worker", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--tasks", default="fold,cofold,subopt")
    parser.add_argument("--warmups", type=int, default=0)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--bootstrap-iterations", type=int, default=5000)
    parser.add_argument("--significance-threshold", type=float, default=0.95)
    parser.add_argument("--expect-avx2", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    requested = [task.strip() for task in args.tasks.split(",") if task.strip()]
    workloads = {
        "fold": fold_workload(),
        "cofold": cofold_workload(),
        "subopt": subopt_workload(),
    }

    new_probe = run_probe(Path(args.new_worker))
    old_probe = run_probe(Path(args.old_worker))
    (out_dir / "new_probe.json").write_text(json.dumps(new_probe, indent=2), encoding="utf-8")
    (out_dir / "old_probe.json").write_text(json.dumps(old_probe, indent=2), encoding="utf-8")

    if args.expect_avx2:
        compiled = set(new_probe.get("compiled_slices", []))
        runtime = set(new_probe.get("runtime_features", []))
        if ("AVX2" not in compiled) or ("AVX2" not in runtime):
            raise SystemExit("AVX2 expected but new worker does not report compiled/runtime AVX2 support")

    overall = {
        "new_probe": new_probe,
        "old_probe": old_probe,
    }

    for task in requested:
        old_rows = run_worker(Path(args.old_worker), task, workloads[task], args.repeats, args.warmups)
        new_rows = run_worker(Path(args.new_worker), task, workloads[task], args.repeats, args.warmups)

        if task == "subopt":
            diffs, rows = compare_subopt(old_rows, new_rows)
        else:
            diffs, rows = compare_fold_like(old_rows, new_rows)

        write_csv(out_dir / f"capi_{task}_compare.csv", rows)
        summary = summarize(task,
                            old_rows,
                            new_rows,
                            diffs,
                            args.bootstrap_iterations,
                            args.significance_threshold)
        (out_dir / f"capi_{task}_compare.summary.json").write_text(json.dumps(summary, indent=2),
                                                                   encoding="utf-8")
        overall[task] = summary

    (out_dir / "capi_compare.summary.json").write_text(json.dumps(overall, indent=2), encoding="utf-8")
    print(json.dumps(overall, indent=2))


if __name__ == "__main__":
    main()

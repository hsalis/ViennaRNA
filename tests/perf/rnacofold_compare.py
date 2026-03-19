#!/usr/bin/env python3

import argparse
import csv
import json
import random
import re
import statistics
import subprocess
import time
from collections import defaultdict
from pathlib import Path


ALPHABET = "ACGU"
DEFAULT_LENGTHS = (50, 100, 200)
DEFAULT_NUM_DIMERS = 30
DEFAULT_REPEATS = 5
DEFAULT_SEED = 1729

RESULT_RE = re.compile(r"^([().&]+)\s+\(\s*([-+]?\d+(?:\.\d+)?)\)")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare RNAcofold predictions and runtime across two ViennaRNA builds."
    )
    parser.add_argument("--new-bin", required=True, help="Path to new RNAcofold binary")
    parser.add_argument("--old-bin", required=True, help="Path to old RNAcofold binary")
    parser.add_argument(
        "--lengths",
        default="50,100,200",
        help="Comma-separated strand lengths to sample from",
    )
    parser.add_argument(
        "--num-dimers",
        type=int,
        default=DEFAULT_NUM_DIMERS,
        help="Number of randomized dimers to compare",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=DEFAULT_REPEATS,
        help="Whole-workload timing repeats per binary",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_SEED,
        help="Random seed",
    )
    parser.add_argument(
        "--output-prefix",
        default="tests/perf/rnacofold_compare",
        help="Output prefix for csv/json/summary files",
    )
    return parser.parse_args()


def random_rna(rng, length):
    return "".join(rng.choice(ALPHABET) for _ in range(length))


def generate_dimers(lengths, num_dimers, seed):
    rng = random.Random(seed)
    combos = [(a, b) for a in lengths for b in lengths]
    dimers = []

    for idx in range(num_dimers):
        len1, len2 = combos[idx % len(combos)]
        dimers.append(
            {
                "id": f"dimer_{idx + 1:02d}",
                "len1": len1,
                "len2": len2,
                "sequence1": random_rna(rng, len1),
                "sequence2": random_rna(rng, len2),
            }
        )

    return dimers


def parse_result(stdout):
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    if len(lines) < 2:
      raise RuntimeError(f"Unexpected RNAcofold output:\n{stdout}")

    match = RESULT_RE.match(lines[-1])
    if not match:
        raise RuntimeError(f"Failed to parse RNAcofold result line: {lines[-1]!r}")

    return match.group(1), float(match.group(2))


def run_case(binary, seq1, seq2):
    payload = f"{seq1}&{seq2}\n"
    start = time.perf_counter_ns()
    proc = subprocess.run(
        [binary, "--noPS"],
        input=payload,
        text=True,
        capture_output=True,
        check=True,
    )
    elapsed = time.perf_counter_ns() - start
    structure, energy = parse_result(proc.stdout)
    return structure, energy, elapsed


def collect_predictions(binary, dimers):
    rows = []
    for dimer in dimers:
        structure, energy, elapsed = run_case(
            binary,
            dimer["sequence1"],
            dimer["sequence2"],
        )
        rows.append(
            {
                **dimer,
                "joined_sequence": f"{dimer['sequence1']}&{dimer['sequence2']}",
                "structure": structure,
                "energy": energy,
                "wall_ns": elapsed,
            }
        )
    return rows


def benchmark_workload(binary, dimers, repeats):
    samples = []

    if dimers:
        run_case(binary, dimers[0]["sequence1"], dimers[0]["sequence2"])

    for _ in range(repeats):
        start = time.perf_counter_ns()
        for dimer in dimers:
            run_case(binary, dimer["sequence1"], dimer["sequence2"])
        samples.append(time.perf_counter_ns() - start)

    return {
        "samples_ns": samples,
        "median_ns": int(statistics.median(samples)),
        "best_ns": min(samples),
        "ns_per_dimer": statistics.median(samples) / len(dimers) if dimers else 0.0,
    }


def compare_rows(new_rows, old_rows):
    old_by_id = {row["id"]: row for row in old_rows}
    comparison = []

    for new_row in new_rows:
        old_row = old_by_id[new_row["id"]]
        delta = new_row["energy"] - old_row["energy"]
        comparison.append(
            {
                "id": new_row["id"],
                "len1": new_row["len1"],
                "len2": new_row["len2"],
                "sequence1": new_row["sequence1"],
                "sequence2": new_row["sequence2"],
                "new_structure": new_row["structure"],
                "new_energy": new_row["energy"],
                "old_structure": old_row["structure"],
                "old_energy": old_row["energy"],
                "same_structure": new_row["structure"] == old_row["structure"],
                "same_energy": abs(delta) < 1e-6,
                "energy_delta": delta,
                "new_wall_ns": new_row["wall_ns"],
                "old_wall_ns": old_row["wall_ns"],
                "new_over_old_runtime": (
                    new_row["wall_ns"] / old_row["wall_ns"] if old_row["wall_ns"] else 0.0
                ),
                "different": (
                    (new_row["structure"] != old_row["structure"]) or (abs(delta) >= 1e-6)
                ),
            }
        )

    return comparison


def summarize(comparison, new_bench, old_bench, seed):
    diffs = [row for row in comparison if row["different"]]
    per_combo = defaultdict(list)

    for row in comparison:
        per_combo[(row["len1"], row["len2"])].append(row)

    combo_summary = {}
    for combo, rows in sorted(per_combo.items()):
        combo_summary[f"{combo[0]}x{combo[1]}"] = {
            "count": len(rows),
            "diff_count": sum(1 for row in rows if row["different"]),
            "median_new_wall_ns": int(statistics.median(row["new_wall_ns"] for row in rows)),
            "median_old_wall_ns": int(statistics.median(row["old_wall_ns"] for row in rows)),
            "median_new_over_old_runtime": statistics.median(
                row["new_over_old_runtime"] for row in rows
            ),
            "max_abs_energy_delta": max(abs(row["energy_delta"]) for row in rows),
        }

    return {
        "seed": seed,
        "total_dimers": len(comparison),
        "total_differences": len(diffs),
        "workload_runtime": {
            "new_median_ns": new_bench["median_ns"],
            "new_best_ns": new_bench["best_ns"],
            "new_ns_per_dimer": new_bench["ns_per_dimer"],
            "old_median_ns": old_bench["median_ns"],
            "old_best_ns": old_bench["best_ns"],
            "old_ns_per_dimer": old_bench["ns_per_dimer"],
            "new_over_old_ratio": (
                new_bench["median_ns"] / old_bench["median_ns"] if old_bench["median_ns"] else 0.0
            ),
        },
        "combo_summary": combo_summary,
        "differences": diffs,
    }


def write_csv(path, rows):
    fieldnames = [
        "id",
        "len1",
        "len2",
        "sequence1",
        "sequence2",
        "new_structure",
        "new_energy",
        "old_structure",
        "old_energy",
        "same_structure",
        "same_energy",
        "energy_delta",
        "new_wall_ns",
        "old_wall_ns",
        "new_over_old_runtime",
        "different",
    ]

    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_summary(path, summary):
    lines = [
        f"seed: {summary['seed']}",
        f"total_dimers: {summary['total_dimers']}",
        f"total_differences: {summary['total_differences']}",
        "",
        "workload_runtime:",
        (
            f"  new_median_ns={summary['workload_runtime']['new_median_ns']} "
            f"old_median_ns={summary['workload_runtime']['old_median_ns']} "
            f"ratio={summary['workload_runtime']['new_over_old_ratio']:.4f}"
        ),
        "",
        "per_length_pair:",
    ]

    for combo, info in summary["combo_summary"].items():
        lines.append(
            f"  {combo}: count={info['count']} diff_count={info['diff_count']} "
            f"median_new_ns={info['median_new_wall_ns']} "
            f"median_old_ns={info['median_old_wall_ns']} "
            f"median_ratio={info['median_new_over_old_runtime']:.4f} "
            f"max_abs_energy_delta={info['max_abs_energy_delta']:.6f}"
        )

    lines.append("")
    if summary["differences"]:
        lines.append("differences:")
        for row in summary["differences"]:
            lines.append(
                f"  {row['id']} {row['len1']}x{row['len2']}: "
                f"new={row['new_structure']} ({row['new_energy']:.2f}) "
                f"old={row['old_structure']} ({row['old_energy']:.2f}) "
                f"delta={row['energy_delta']:.6f}"
            )
    else:
        lines.append("differences: none")

    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    args = parse_args()
    lengths = [int(value.strip()) for value in args.lengths.split(",") if value.strip()]
    output_prefix = Path(args.output_prefix)
    output_prefix.parent.mkdir(parents=True, exist_ok=True)

    dimers = generate_dimers(lengths, args.num_dimers, args.seed)
    new_rows = collect_predictions(args.new_bin, dimers)
    old_rows = collect_predictions(args.old_bin, dimers)
    comparison = compare_rows(new_rows, old_rows)

    new_bench = benchmark_workload(args.new_bin, dimers, args.repeats)
    old_bench = benchmark_workload(args.old_bin, dimers, args.repeats)
    summary = summarize(comparison, new_bench, old_bench, args.seed)

    csv_path = output_prefix.with_suffix(".csv")
    json_path = output_prefix.with_suffix(".json")
    summary_path = output_prefix.with_suffix(".summary.txt")

    write_csv(csv_path, comparison)
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_summary(summary_path, summary)

    print(f"csv={csv_path}")
    print(f"json={json_path}")
    print(f"summary={summary_path}")
    print(f"total_dimers={summary['total_dimers']}")
    print(f"total_differences={summary['total_differences']}")
    print(
        "runtime_new_median_ns="
        f"{summary['workload_runtime']['new_median_ns']}"
    )
    print(
        "runtime_old_median_ns="
        f"{summary['workload_runtime']['old_median_ns']}"
    )
    print(
        "runtime_new_over_old_ratio="
        f"{summary['workload_runtime']['new_over_old_ratio']:.4f}"
    )


if __name__ == "__main__":
    main()

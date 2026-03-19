#!/usr/bin/env python3

import argparse
import csv
import json
import os
import random
import statistics
import subprocess
import sys
import tempfile
from pathlib import Path


ALPHABET = "ACGU"


CHILD_CODE = r"""
import json
import sys

import RNA

input_path = sys.argv[1]
output_path = sys.argv[2]

with open(input_path, "r", encoding="utf-8") as handle:
    payload = json.load(handle)

results = []
for entry in payload["sequences"]:
    structure, energy = RNA.fold(entry["sequence"])
    results.append({
        "id": entry["id"],
        "length": entry["length"],
        "sequence": entry["sequence"],
        "structure": structure,
        "energy": float(energy),
    })

with open(output_path, "w", encoding="utf-8") as handle:
    json.dump({
        "version": getattr(RNA, "__version__", "unknown"),
        "results": results,
    }, handle)
"""


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare RNA.fold predictions across two ViennaRNA Python builds."
    )
    parser.add_argument("--new-module-dir", required=True, help="Path to new RNA Python module dir")
    parser.add_argument("--old-module-dir", required=True, help="Path to old RNA Python module dir")
    parser.add_argument(
        "--lengths",
        default="100,500,1000",
        help="Comma-separated sequence lengths",
    )
    parser.add_argument(
        "--seqs-per-length",
        type=int,
        default=30,
        help="Number of randomized sequences per length",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1729,
        help="Random seed",
    )
    parser.add_argument(
        "--output-prefix",
        default="tests/perf/python_rnafold_compare",
        help="Output prefix for csv/json/summary files",
    )
    return parser.parse_args()


def generate_sequences(lengths, seqs_per_length, seed):
    rng = random.Random(seed)
    sequences = []
    for length in lengths:
        for index in range(seqs_per_length):
            seq = "".join(rng.choice(ALPHABET) for _ in range(length))
            sequences.append(
                {
                    "id": f"len{length}_{index + 1:02d}",
                    "length": length,
                    "sequence": seq,
                }
            )
    return sequences


def run_module(module_dir, sequences):
    payload = {"sequences": sequences}
    env = os.environ.copy()
    env["PYTHONPATH"] = (
        str(Path(module_dir).resolve())
        + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
    )

    with tempfile.TemporaryDirectory(prefix="vrna-compare-") as tmpdir:
        input_path = Path(tmpdir) / "input.json"
        output_path = Path(tmpdir) / "output.json"
        input_path.write_text(json.dumps(payload), encoding="utf-8")
        subprocess.run(
            [sys.executable, "-c", CHILD_CODE, str(input_path), str(output_path)],
            check=True,
            env=env,
        )
        return json.loads(output_path.read_text(encoding="utf-8"))


def compare_results(new_payload, old_payload):
    old_by_id = {entry["id"]: entry for entry in old_payload["results"]}
    rows = []
    for new_entry in new_payload["results"]:
        old_entry = old_by_id[new_entry["id"]]
        same_structure = new_entry["structure"] == old_entry["structure"]
        energy_delta = new_entry["energy"] - old_entry["energy"]
        same_energy = abs(energy_delta) < 1e-6
        rows.append(
            {
                "id": new_entry["id"],
                "length": new_entry["length"],
                "sequence": new_entry["sequence"],
                "new_structure": new_entry["structure"],
                "new_energy": new_entry["energy"],
                "old_structure": old_entry["structure"],
                "old_energy": old_entry["energy"],
                "same_structure": same_structure,
                "same_energy": same_energy,
                "energy_delta": energy_delta,
                "different": (not same_structure) or (not same_energy),
            }
        )
    return rows


def write_csv(path, rows):
    fieldnames = [
        "id",
        "length",
        "sequence",
        "new_structure",
        "new_energy",
        "old_structure",
        "old_energy",
        "same_structure",
        "same_energy",
        "energy_delta",
        "different",
    ]
    with open(path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def summarize(rows, new_version, old_version, seed, lengths, seqs_per_length):
    diffs = [row for row in rows if row["different"]]
    by_length = {}
    for length in sorted({row["length"] for row in rows}):
        group = [row for row in rows if row["length"] == length]
        group_diffs = [row for row in group if row["different"]]
        energy_deltas = [abs(row["energy_delta"]) for row in group]
        by_length[length] = {
            "count": len(group),
            "diff_count": len(group_diffs),
            "max_abs_energy_delta": max(energy_deltas) if energy_deltas else 0.0,
            "median_abs_energy_delta": statistics.median(energy_deltas) if energy_deltas else 0.0,
        }

    return {
        "new_version": new_version,
        "old_version": old_version,
        "seed": seed,
        "lengths": lengths,
        "seqs_per_length": seqs_per_length,
        "total_sequences": len(rows),
        "total_differences": len(diffs),
        "length_summary": by_length,
        "differences": diffs,
    }


def write_summary(path, summary):
    lines = []
    lines.append(f"new_version: {summary['new_version']}")
    lines.append(f"old_version: {summary['old_version']}")
    lines.append(f"seed: {summary['seed']}")
    lines.append(f"lengths: {', '.join(str(v) for v in summary['lengths'])}")
    lines.append(f"seqs_per_length: {summary['seqs_per_length']}")
    lines.append(f"total_sequences: {summary['total_sequences']}")
    lines.append(f"total_differences: {summary['total_differences']}")
    lines.append("")
    lines.append("per_length:")
    for length, info in summary["length_summary"].items():
        lines.append(
            f"  {length}: count={info['count']} diff_count={info['diff_count']} "
            f"max_abs_energy_delta={info['max_abs_energy_delta']:.6f} "
            f"median_abs_energy_delta={info['median_abs_energy_delta']:.6f}"
        )
    lines.append("")
    if summary["differences"]:
        lines.append("differences:")
        for row in summary["differences"]:
            lines.append(
                f"  {row['id']}: len={row['length']} "
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

    sequences = generate_sequences(lengths, args.seqs_per_length, args.seed)
    new_payload = run_module(args.new_module_dir, sequences)
    old_payload = run_module(args.old_module_dir, sequences)
    rows = compare_results(new_payload, old_payload)

    csv_path = output_prefix.with_suffix(".csv")
    json_path = output_prefix.with_suffix(".json")
    summary_path = output_prefix.with_suffix(".summary.txt")

    write_csv(csv_path, rows)
    summary = summarize(
        rows,
        new_payload["version"],
        old_payload["version"],
        args.seed,
        lengths,
        args.seqs_per_length,
    )
    write_summary(summary_path, summary)
    json_path.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"csv={csv_path}")
    print(f"json={json_path}")
    print(f"summary={summary_path}")
    print(f"total_sequences={summary['total_sequences']}")
    print(f"total_differences={summary['total_differences']}")
    for length, info in summary["length_summary"].items():
        print(
            f"length={length} count={info['count']} diff_count={info['diff_count']} "
            f"max_abs_energy_delta={info['max_abs_energy_delta']:.6f}"
        )


if __name__ == "__main__":
    main()

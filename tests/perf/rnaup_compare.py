#!/usr/bin/env python3

import argparse
import csv
import json
import random
import statistics
import subprocess
import time
from pathlib import Path


ALPHABET = "ACGU"
LENGTHS = (50, 100, 200)


def synthetic_sequence(length: int, seed: int) -> str:
    rng = random.Random(seed ^ length)
    return "".join(rng.choice(ALPHABET) for _ in range(length))


def run_rnaup(binary: Path, seq1: str, seq2: str) -> str:
    proc = subprocess.run(
        [str(binary), "-b", "-o", "--no_header"],
        input=f"{seq1}\n{seq2}\n",
        text=True,
        capture_output=True,
        check=True,
    )
    return proc.stdout.strip()


def time_rnaup(binary: Path, seq1: str, seq2: str, repeats: int) -> int:
    samples = []
    for _ in range(repeats):
      start = time.perf_counter_ns()
      run_rnaup(binary, seq1, seq2)
      samples.append(time.perf_counter_ns() - start)
    return int(statistics.median(samples))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--new-binary", required=True)
    parser.add_argument("--old-binary", required=True)
    parser.add_argument("--cases", type=int, default=30)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--seed", type=int, default=1729)
    parser.add_argument("--csv-out", required=True)
    parser.add_argument("--json-out", required=True)
    args = parser.parse_args()

    new_binary = Path(args.new_binary)
    old_binary = Path(args.old_binary)
    csv_out = Path(args.csv_out)
    json_out = Path(args.json_out)

    rng = random.Random(args.seed)
    rows = []

    for case_id in range(1, args.cases + 1):
        len1 = rng.choice(LENGTHS)
        len2 = rng.choice(LENGTHS)
        seq1 = synthetic_sequence(len1, args.seed + case_id * 11)
        seq2 = synthetic_sequence(len2, args.seed + case_id * 17)

        new_out = run_rnaup(new_binary, seq1, seq2)
        old_out = run_rnaup(old_binary, seq1, seq2)
        new_ns = time_rnaup(new_binary, seq1, seq2, args.repeats)
        old_ns = time_rnaup(old_binary, seq1, seq2, args.repeats)

        rows.append(
            {
                "case_id": f"rnaup_{case_id:03d}",
                "len1": len1,
                "len2": len2,
                "seq1": seq1,
                "seq2": seq2,
                "new_output": new_out,
                "old_output": old_out,
                "match": int(new_out == old_out),
                "new_wall_ns": new_ns,
                "old_wall_ns": old_ns,
                "runtime_ratio_new_over_old": (float(new_ns) / float(old_ns)) if old_ns else 0.0,
            }
        )

    with csv_out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    total_new = sum(row["new_wall_ns"] for row in rows)
    total_old = sum(row["old_wall_ns"] for row in rows)
    summary = {
        "cases": len(rows),
        "differences": sum(1 for row in rows if not row["match"]),
        "matches": sum(1 for row in rows if row["match"]),
        "total_new_wall_ns": total_new,
        "total_old_wall_ns": total_old,
        "runtime_ratio_new_over_old": (float(total_new) / float(total_old)) if total_old else 0.0,
    }

    with json_out.open("w") as f:
        json.dump(summary, f, indent=2, sort_keys=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3

import argparse
import csv
import json
import pathlib
import random
import re
import statistics
import subprocess
import tempfile
import time


TERMINAL_RE = re.compile(r"(-?\d+(?:\.\d+)?)\s+(O|X\d+)\s*$")
TRAJECTORY_RE = re.compile(
    r"^([().]+)\s+(-?\d+(?:\.\d+)?)\s+(-?\d+(?:\.\d+)?)\s*(O|X\d+)?\s*$"
)
RNA_ALPHABET = "ACGU"


def make_seed_triplet(rng):
    return tuple(rng.randint(1, 65535) for _ in range(3))


def random_sequence(rng, length):
    return "".join(rng.choice(RNA_ALPHABET) for _ in range(length))


def write_input_file(sequence):
    tmp = tempfile.NamedTemporaryFile("w", prefix="kinfold-random-", suffix=".in", delete=False)
    tmp.write(sequence)
    tmp.write("\n")
    tmp.close()
    return tmp.name


def parse_trajectory_records(stdout):
    records = []
    terminal_events = []

    for raw_line in stdout.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue

        match = TRAJECTORY_RE.match(line)
        if not match:
            continue

        marker = match.group(4)
        record = {
            "line": line,
            "structure": match.group(1),
            "energy": float(match.group(2)),
            "time": float(match.group(3)),
            "marker": marker,
        }
        records.append(record)

        if marker:
            terminal_events.append(
                {
                    "line": line,
                    "time": record["time"],
                    "marker": marker,
                }
            )

    return records, terminal_events


def canonical_most_stable(records):
    if not records:
        return None

    min_energy = min(record["energy"] for record in records)
    stable_structures = sorted(
        {
            record["structure"]
            for record in records
            if record["energy"] == min_energy
        }
    )
    structure = stable_structures[0]
    occurrences = sum(
        1
        for record in records
        if record["energy"] == min_energy and record["structure"] == structure
    )
    return {
        "structure": structure,
        "energy": min_energy,
        "degeneracy": len(stable_structures),
        "occurrences": occurrences,
        "all_structures": stable_structures,
    }


def run_case(binary, input_file, num_trajectories, seed_triplet, extra_args):
    cmd = [
        binary,
        "--fpt",
        "--num",
        str(num_trajectories),
        "--seed",
        f"{seed_triplet[0]}={seed_triplet[1]}={seed_triplet[2]}",
    ] + extra_args

    with tempfile.TemporaryDirectory(prefix="kinfold-random-run-") as tmpdir:
        start = time.perf_counter_ns()
        with open(input_file, "rb") as handle:
            proc = subprocess.run(
                cmd,
                stdin=handle,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
                cwd=tmpdir,
            )
        wall_ns = time.perf_counter_ns() - start

    stdout = proc.stdout.decode("utf-8")
    records, terminal_events = parse_trajectory_records(stdout)

    if not terminal_events:
        raise RuntimeError(f"{binary} produced no terminal events to parse")

    fold_times = [event["time"] for event in terminal_events]
    markers = [event["marker"] for event in terminal_events]
    most_stable = canonical_most_stable(records)

    return {
        "wall_ns": wall_ns,
        "events": terminal_events,
        "fold_times": fold_times,
        "markers": markers,
        "terminal_count": len(terminal_events),
        "success_count": sum(1 for marker in markers if marker.startswith("X")),
        "timeout_count": sum(1 for marker in markers if marker == "O"),
        "median_fold_time": statistics.median(fold_times),
        "mean_fold_time": statistics.fmean(fold_times),
        "min_fold_time": min(fold_times),
        "max_fold_time": max(fold_times),
        "most_stable": most_stable,
    }


def summarize_single_cases(cases):
    by_length = {}
    for case in cases:
        bucket = by_length.setdefault(
            case["length"],
            {
                "length": case["length"],
                "case_count": 0,
                "trajectory_count": 0,
                "success_count": 0,
                "timeout_count": 0,
                "wall_ns": [],
                "mean_fold_time": [],
            },
        )
        bucket["case_count"] += 1
        bucket["trajectory_count"] += case["terminal_count"]
        bucket["success_count"] += case["success_count"]
        bucket["timeout_count"] += case["timeout_count"]
        bucket["wall_ns"].append(case["wall_ns"])
        bucket["mean_fold_time"].append(case["mean_fold_time"])

    summary_by_length = []
    for length in sorted(by_length):
        bucket = by_length[length]
        summary_by_length.append(
            {
                "length": length,
                "case_count": bucket["case_count"],
                "trajectory_count": bucket["trajectory_count"],
                "success_count": bucket["success_count"],
                "timeout_count": bucket["timeout_count"],
                "success_rate": (
                    float(bucket["success_count"]) / float(bucket["trajectory_count"])
                    if bucket["trajectory_count"]
                    else None
                ),
                "median_wall_ns": statistics.median(bucket["wall_ns"]),
                "mean_wall_ns": statistics.fmean(bucket["wall_ns"]),
                "mean_of_case_mean_fold_time": statistics.fmean(bucket["mean_fold_time"]),
            }
        )

    all_wall_ns = [case["wall_ns"] for case in cases]
    all_fold_times = [time_value for case in cases for time_value in case["fold_times"]]
    total_trajectories = sum(case["terminal_count"] for case in cases)
    total_success = sum(case["success_count"] for case in cases)
    total_timeout = sum(case["timeout_count"] for case in cases)

    overall = {
        "case_count": len(cases),
        "trajectory_count": total_trajectories,
        "success_count": total_success,
        "timeout_count": total_timeout,
        "success_rate": (float(total_success) / float(total_trajectories)) if total_trajectories else None,
        "median_wall_ns": statistics.median(all_wall_ns) if all_wall_ns else None,
        "mean_wall_ns": statistics.fmean(all_wall_ns) if all_wall_ns else None,
        "mean_of_case_mean_fold_time": statistics.fmean(case["mean_fold_time"] for case in cases)
        if cases
        else None,
        "mean_fold_time": statistics.fmean(all_fold_times) if all_fold_times else None,
        "min_fold_time": min(all_fold_times) if all_fold_times else None,
        "max_fold_time": max(all_fold_times) if all_fold_times else None,
    }

    return summary_by_length, overall


def summarize_compare_cases(cases):
    by_length = {}
    for case in cases:
        bucket = by_length.setdefault(
            case["length"],
            {
                "length": case["length"],
                "case_count": 0,
                "old_wall_ns": [],
                "new_wall_ns": [],
                "old_mean_fold_time": [],
                "new_mean_fold_time": [],
                "stable_structure_match_count": 0,
                "stable_energy_match_count": 0,
            },
        )
        bucket["case_count"] += 1
        bucket["old_wall_ns"].append(case["old"]["wall_ns"])
        bucket["new_wall_ns"].append(case["new"]["wall_ns"])
        bucket["old_mean_fold_time"].append(case["old"]["mean_fold_time"])
        bucket["new_mean_fold_time"].append(case["new"]["mean_fold_time"])
        bucket["stable_structure_match_count"] += int(case["most_stable_structure_match"])
        bucket["stable_energy_match_count"] += int(case["most_stable_energy_match"])

    summary_by_length = []
    for length in sorted(by_length):
        bucket = by_length[length]
        old_total_wall_ns = sum(bucket["old_wall_ns"])
        new_total_wall_ns = sum(bucket["new_wall_ns"])
        summary_by_length.append(
            {
                "length": length,
                "case_count": bucket["case_count"],
                "old_mean_of_case_mean_fold_time": statistics.fmean(bucket["old_mean_fold_time"]),
                "new_mean_of_case_mean_fold_time": statistics.fmean(bucket["new_mean_fold_time"]),
                "old_mean_wall_ns": statistics.fmean(bucket["old_wall_ns"]),
                "new_mean_wall_ns": statistics.fmean(bucket["new_wall_ns"]),
                "old_total_wall_ns": old_total_wall_ns,
                "new_total_wall_ns": new_total_wall_ns,
                "speedup_old_over_new": (
                    float(old_total_wall_ns) / float(new_total_wall_ns)
                    if new_total_wall_ns
                    else None
                ),
                "most_stable_structure_match_rate": (
                    float(bucket["stable_structure_match_count"]) / float(bucket["case_count"])
                    if bucket["case_count"]
                    else None
                ),
                "most_stable_energy_match_rate": (
                    float(bucket["stable_energy_match_count"]) / float(bucket["case_count"])
                    if bucket["case_count"]
                    else None
                ),
            }
        )

    old_total_wall_ns = sum(case["old"]["wall_ns"] for case in cases)
    new_total_wall_ns = sum(case["new"]["wall_ns"] for case in cases)
    overall = {
        "case_count": len(cases),
        "old_mean_of_case_mean_fold_time": statistics.fmean(case["old"]["mean_fold_time"] for case in cases)
        if cases
        else None,
        "new_mean_of_case_mean_fold_time": statistics.fmean(case["new"]["mean_fold_time"] for case in cases)
        if cases
        else None,
        "old_mean_wall_ns": statistics.fmean(case["old"]["wall_ns"] for case in cases) if cases else None,
        "new_mean_wall_ns": statistics.fmean(case["new"]["wall_ns"] for case in cases) if cases else None,
        "old_total_wall_ns": old_total_wall_ns,
        "new_total_wall_ns": new_total_wall_ns,
        "speedup_old_over_new": (
            float(old_total_wall_ns) / float(new_total_wall_ns) if new_total_wall_ns else None
        ),
        "most_stable_structure_match_rate": (
            statistics.fmean(int(case["most_stable_structure_match"]) for case in cases) if cases else None
        ),
        "most_stable_energy_match_rate": (
            statistics.fmean(int(case["most_stable_energy_match"]) for case in cases) if cases else None
        ),
    }

    return summary_by_length, overall


def write_single_csv(path, cases):
    fieldnames = [
        "case_id",
        "length",
        "sequence",
        "seed_triplet",
        "terminal_count",
        "success_count",
        "timeout_count",
        "success_rate",
        "wall_ns",
        "mean_fold_time",
        "median_fold_time",
        "min_fold_time",
        "max_fold_time",
        "most_stable_energy",
        "most_stable_structure",
    ]

    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for case in cases:
            writer.writerow(
                {
                    "case_id": case["case_id"],
                    "length": case["length"],
                    "sequence": case["sequence"],
                    "seed_triplet": "=".join(str(x) for x in case["seed_triplet"]),
                    "terminal_count": case["terminal_count"],
                    "success_count": case["success_count"],
                    "timeout_count": case["timeout_count"],
                    "success_rate": (
                        float(case["success_count"]) / float(case["terminal_count"])
                        if case["terminal_count"]
                        else None
                    ),
                    "wall_ns": case["wall_ns"],
                    "mean_fold_time": case["mean_fold_time"],
                    "median_fold_time": case["median_fold_time"],
                    "min_fold_time": case["min_fold_time"],
                    "max_fold_time": case["max_fold_time"],
                    "most_stable_energy": case["most_stable"]["energy"],
                    "most_stable_structure": case["most_stable"]["structure"],
                }
            )


def write_compare_csv(path, cases):
    fieldnames = [
        "case_id",
        "length",
        "sequence",
        "seed_triplet",
        "old_mean_fold_time",
        "new_mean_fold_time",
        "old_wall_ns",
        "new_wall_ns",
        "speedup_old_over_new",
        "old_most_stable_energy",
        "new_most_stable_energy",
        "old_most_stable_structure",
        "new_most_stable_structure",
        "most_stable_energy_match",
        "most_stable_structure_match",
    ]

    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for case in cases:
            writer.writerow(
                {
                    "case_id": case["case_id"],
                    "length": case["length"],
                    "sequence": case["sequence"],
                    "seed_triplet": "=".join(str(x) for x in case["seed_triplet"]),
                    "old_mean_fold_time": case["old"]["mean_fold_time"],
                    "new_mean_fold_time": case["new"]["mean_fold_time"],
                    "old_wall_ns": case["old"]["wall_ns"],
                    "new_wall_ns": case["new"]["wall_ns"],
                    "speedup_old_over_new": case["speedup_old_over_new"],
                    "old_most_stable_energy": case["old"]["most_stable"]["energy"],
                    "new_most_stable_energy": case["new"]["most_stable"]["energy"],
                    "old_most_stable_structure": case["old"]["most_stable"]["structure"],
                    "new_most_stable_structure": case["new"]["most_stable"]["structure"],
                    "most_stable_energy_match": case["most_stable_energy_match"],
                    "most_stable_structure_match": case["most_stable_structure_match"],
                }
            )


def build_cases(args, extra_args, compare_mode):
    rng = random.Random(args.seed)
    cases = []
    case_index = 0

    for length in range(args.min_length, args.max_length + 1):
        for sample_index in range(args.samples_per_length):
            sequence = random_sequence(rng, length)
            seed_triplet = make_seed_triplet(rng)
            input_path = write_input_file(sequence)

            try:
                if compare_mode:
                    for _ in range(args.warmups):
                        run_case(args.old_binary, input_path, args.num_trajectories, seed_triplet, extra_args)
                        run_case(args.new_binary, input_path, args.num_trajectories, seed_triplet, extra_args)

                    old_result = run_case(
                        args.old_binary,
                        input_path,
                        args.num_trajectories,
                        seed_triplet,
                        extra_args,
                    )
                    new_result = run_case(
                        args.new_binary,
                        input_path,
                        args.num_trajectories,
                        seed_triplet,
                        extra_args,
                    )
                    case_payload = {
                        "case_id": f"L{length:03d}_S{sample_index:03d}",
                        "length": length,
                        "sample_index": sample_index,
                        "sequence": sequence,
                        "seed_triplet": seed_triplet,
                        "old": old_result,
                        "new": new_result,
                        "speedup_old_over_new": (
                            float(old_result["wall_ns"]) / float(new_result["wall_ns"])
                            if new_result["wall_ns"]
                            else None
                        ),
                        "most_stable_energy_match": (
                            old_result["most_stable"]["energy"] == new_result["most_stable"]["energy"]
                        ),
                        "most_stable_structure_match": (
                            old_result["most_stable"]["structure"]
                            == new_result["most_stable"]["structure"]
                        ),
                    }
                else:
                    for _ in range(args.warmups):
                        run_case(args.binary, input_path, args.num_trajectories, seed_triplet, extra_args)

                    result = run_case(
                        args.binary,
                        input_path,
                        args.num_trajectories,
                        seed_triplet,
                        extra_args,
                    )
                    case_payload = {
                        "case_id": f"L{length:03d}_S{sample_index:03d}",
                        "length": length,
                        "sample_index": sample_index,
                        "sequence": sequence,
                        "seed_triplet": seed_triplet,
                        **result,
                    }
            finally:
                pathlib.Path(input_path).unlink(missing_ok=True)

            cases.append(case_payload)
            case_index += 1

    return cases, case_index


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary")
    parser.add_argument("--old-binary")
    parser.add_argument("--new-binary")
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--output-csv")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--min-length", type=int, default=20)
    parser.add_argument("--max-length", type=int, default=100)
    parser.add_argument("--samples-per-length", type=int, default=3)
    parser.add_argument("--num-trajectories", type=int, default=1000)
    parser.add_argument("--warmups", type=int, default=0)
    args, extra_args = parser.parse_known_args()

    compare_mode = bool(args.old_binary or args.new_binary)
    if compare_mode:
        if args.binary:
            raise ValueError("use either --binary or both --old-binary/--new-binary")
        if not (args.old_binary and args.new_binary):
            raise ValueError("compare mode requires both --old-binary and --new-binary")
    else:
        if not args.binary:
            raise ValueError("single-binary mode requires --binary")

    if args.min_length < 1:
        raise ValueError("--min-length must be positive")
    if args.max_length < args.min_length:
        raise ValueError("--max-length must be >= --min-length")
    if args.samples_per_length < 1:
        raise ValueError("--samples-per-length must be >= 1")
    if args.num_trajectories < 1:
        raise ValueError("--num-trajectories must be >= 1")

    cases, case_count = build_cases(args, extra_args, compare_mode)

    if compare_mode:
        summary_by_length, overall = summarize_compare_cases(cases)
        payload = {
            "mode": "compare",
            "old_binary": args.old_binary,
            "new_binary": args.new_binary,
            "seed": args.seed,
            "min_length": args.min_length,
            "max_length": args.max_length,
            "samples_per_length": args.samples_per_length,
            "num_trajectories": args.num_trajectories,
            "warmups": args.warmups,
            "extra_args": extra_args,
            "case_count": case_count,
            "cases": cases,
            "summary_by_length": summary_by_length,
            "overall": overall,
        }
    else:
        summary_by_length, overall = summarize_single_cases(cases)
        payload = {
            "mode": "single",
            "binary": args.binary,
            "seed": args.seed,
            "min_length": args.min_length,
            "max_length": args.max_length,
            "samples_per_length": args.samples_per_length,
            "num_trajectories": args.num_trajectories,
            "warmups": args.warmups,
            "extra_args": extra_args,
            "case_count": case_count,
            "cases": cases,
            "summary_by_length": summary_by_length,
            "overall": overall,
        }

    pathlib.Path(args.output_json).write_text(json.dumps(payload, indent=2, sort_keys=True))

    if args.output_csv:
        if compare_mode:
            write_compare_csv(args.output_csv, cases)
        else:
            write_single_csv(args.output_csv, cases)


if __name__ == "__main__":
    main()

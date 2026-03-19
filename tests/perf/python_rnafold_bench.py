#!/usr/bin/env python3

import argparse
import gc
import importlib
import random
import statistics
import sys
import time
from pathlib import Path


DEFAULT_LENGTHS = (100, 500, 1000)
DEFAULT_SEQS_PER_LENGTH = 50
DEFAULT_REPEATS = 5
ALPHABET = "ACGU"


def parse_args():
    parser = argparse.ArgumentParser(
        description="Benchmark Python RNAfold bindings on randomized sequences."
    )
    parser.add_argument(
        "--module-dir",
        type=Path,
        default=None,
        help="Directory containing the built Python RNA package",
    )
    parser.add_argument(
        "--lengths",
        type=int,
        nargs="+",
        default=list(DEFAULT_LENGTHS),
        help="Sequence lengths to benchmark",
    )
    parser.add_argument(
        "--seqs-per-length",
        type=int,
        default=DEFAULT_SEQS_PER_LENGTH,
        help="Number of randomized sequences per length",
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=DEFAULT_REPEATS,
        help="Timing repeats per method and length",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=1,
        help="Random seed",
    )
    parser.add_argument(
        "--csv",
        action="store_true",
        help="Print CSV only",
    )
    return parser.parse_args()


def candidate_module_dirs(root: Path):
  return [
    root / "build-codex-nompfr" / "interfaces" / "Python",
    root / "build-codex-validate" / "interfaces" / "Python",
    root / "interfaces" / "Python",
  ]


def import_rna(module_dir: Path | None):
    root = Path(__file__).resolve().parents[2]
    candidates = [module_dir] if module_dir else candidate_module_dirs(root)

    for candidate in candidates:
        if not candidate:
            continue

        candidate = candidate.resolve()
        if not candidate.exists():
            continue

        sys.path.insert(0, str(candidate))
        try:
            return importlib.import_module("RNA"), candidate
        except ImportError:
            sys.path.pop(0)

    raise ImportError("Could not import local RNA Python module")


def random_sequences(rng: random.Random, length: int, count: int):
    return [
        "".join(rng.choice(ALPHABET) for _ in range(length))
        for _ in range(count)
    ]


def verify_results(rna_module, sequences):
    for seq in sequences:
        structure_a, energy_a = rna_module.fold(seq)
        fc = rna_module.fold_compound(seq)
        structure_b, energy_b = fc.mfe()

        if structure_a != structure_b:
            raise RuntimeError(f"Structure mismatch for sequence length {len(seq)}")

        if abs(energy_a - energy_b) > 1e-5:
            raise RuntimeError(f"Energy mismatch for sequence length {len(seq)}")


def bench_fold(rna_module, sequences):
    start = time.perf_counter_ns()
    for seq in sequences:
        rna_module.fold(seq)
    return time.perf_counter_ns() - start


def bench_fold_compound_mfe(rna_module, sequences):
    start = time.perf_counter_ns()
    for seq in sequences:
        fc = rna_module.fold_compound(seq)
        fc.mfe()
    return time.perf_counter_ns() - start


def run_method(name, func, rna_module, sequences, repeats):
    samples = []

    for _ in range(repeats):
        gc.collect()
        gc.disable()
        try:
            samples.append(func(rna_module, sequences))
        finally:
            gc.enable()

    median_ns = int(statistics.median(samples))
    best_ns = min(samples)

    return {
        "method": name,
        "samples": samples,
        "median_ns": median_ns,
        "best_ns": best_ns,
        "ns_per_sequence": median_ns / len(sequences),
    }


def main():
    args = parse_args()
    repo_root = Path(__file__).resolve().parents[2]
    rna_module, module_path = import_rna(args.module_dir)
    rng = random.Random(args.seed)

    rows = []

    for length in args.lengths:
        sequences = random_sequences(rng, length, args.seqs_per_length)
        verify_results(rna_module, sequences[: min(3, len(sequences))])

        # warm up both call paths once before timing
        bench_fold(rna_module, sequences[:1])
        bench_fold_compound_mfe(rna_module, sequences[:1])

        fold_stats = run_method("RNA.fold", bench_fold, rna_module, sequences, args.repeats)
        fc_stats = run_method(
            "fold_compound(seq).mfe()",
            bench_fold_compound_mfe,
            rna_module,
            sequences,
            args.repeats,
        )

        ratio = fc_stats["median_ns"] / fold_stats["median_ns"] if fold_stats["median_ns"] else 0.0

        rows.append(
            {
                "length": length,
                "count": len(sequences),
                "fold_median_ns": fold_stats["median_ns"],
                "fold_best_ns": fold_stats["best_ns"],
                "fold_ns_per_seq": fold_stats["ns_per_sequence"],
                "fc_median_ns": fc_stats["median_ns"],
                "fc_best_ns": fc_stats["best_ns"],
                "fc_ns_per_seq": fc_stats["ns_per_sequence"],
                "fc_over_fold_ratio": ratio,
            }
        )

    if args.csv:
        print(f"# module_dir={module_path}")
        print("length,count,fold_median_ns,fold_best_ns,fold_ns_per_seq,fc_median_ns,fc_best_ns,fc_ns_per_seq,fc_over_fold_ratio")
        for row in rows:
            print(
                f"{row['length']},{row['count']},"
                f"{row['fold_median_ns']},{row['fold_best_ns']},{row['fold_ns_per_seq']:.2f},"
                f"{row['fc_median_ns']},{row['fc_best_ns']},{row['fc_ns_per_seq']:.2f},"
                f"{row['fc_over_fold_ratio']:.4f}"
            )
        return

    print(f"module_dir: {module_path}")
    print("length count RNA.fold median_ns fold_compound.mfe median_ns ratio")
    for row in rows:
        print(
            f"{row['length']:>6} {row['count']:>5} "
            f"{row['fold_median_ns']:>16} "
            f"{row['fc_median_ns']:>28} "
            f"{row['fc_over_fold_ratio']:>7.4f}"
        )


if __name__ == "__main__":
    main()

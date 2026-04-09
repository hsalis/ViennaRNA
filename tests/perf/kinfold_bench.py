#!/usr/bin/env python3

import argparse
import json
import pathlib
import statistics
import subprocess
import tempfile
import time


def run_once(binary, input_file, extra_args):
    with tempfile.TemporaryDirectory(prefix="kinfold-bench-") as tmpdir:
        cmd = [binary] + extra_args
        start = time.perf_counter_ns()
        subprocess.run(
            cmd,
            stdin=open(input_file, "rb"),
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=True,
            cwd=tmpdir,
        )
        return time.perf_counter_ns() - start


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True)
    parser.add_argument("--input-file", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--warmups", type=int, default=1)
    args, extra_args = parser.parse_known_args()

    for _ in range(args.warmups):
        run_once(args.binary, args.input_file, extra_args)

    timings = [run_once(args.binary, args.input_file, extra_args) for _ in range(args.repeats)]
    payload = {
        "binary": args.binary,
        "input_file": args.input_file,
        "extra_args": extra_args,
        "repeats": args.repeats,
        "warmups": args.warmups,
        "timings_ns": timings,
        "median_ns": statistics.median(timings),
        "total_ns": sum(timings),
    }
    pathlib.Path(args.output).write_text(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

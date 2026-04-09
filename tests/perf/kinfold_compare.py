#!/usr/bin/env python3

import argparse
import json
import os
import pathlib
import re
import statistics
import subprocess
import tempfile
import time

TERMINAL_RE = re.compile(r"(?:^|\s)(O|X\d+)\s*$")


def parse_stdout(stdout):
    lines = [line.strip() for line in stdout.splitlines() if line.strip()]
    markers = []
    for line in lines:
        match = TERMINAL_RE.search(line)
        if match:
            markers.append(match.group(1))
    return {
        "lines": lines,
        "markers": markers,
    }


def run_kinfold(binary, input_path, extra_args, repeats, warmups):
    timings = []
    stdout = ""
    stderr = ""
    log_text = None

    for _ in range(warmups):
      with tempfile.TemporaryDirectory(prefix="kinfold-compare-") as tmpdir:
        log_base = os.path.join(tmpdir, "kinout")
        cmd = [binary, "--log", log_base] + extra_args
        subprocess.run(
            cmd,
            stdin=open(input_path, "rb"),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            cwd=tmpdir,
        )

    for _ in range(repeats):
      with tempfile.TemporaryDirectory(prefix="kinfold-compare-") as tmpdir:
        log_base = os.path.join(tmpdir, "kinout")
        cmd = [binary, "--log", log_base] + extra_args
        start = time.perf_counter_ns()
        proc = subprocess.run(
            cmd,
            stdin=open(input_path, "rb"),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            cwd=tmpdir,
        )
        timings.append(time.perf_counter_ns() - start)
        stdout = proc.stdout.decode("utf-8")
        stderr = proc.stderr.decode("utf-8")
        log_path = log_base + ".log"
        if os.path.exists(log_path):
            log_text = pathlib.Path(log_path).read_text()

    parsed = parse_stdout(stdout)
    return {
        "stdout": stdout,
        "stderr": stderr,
        "log": log_text,
        "parsed": parsed,
        "timings_ns": timings,
        "median_ns": statistics.median(timings),
        "total_ns": sum(timings),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--new-bin", required=True)
    parser.add_argument("--old-bin", required=True)
    parser.add_argument("--input-file", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--warmups", type=int, default=1)
    args, extra_args = parser.parse_known_args()

    old_run = run_kinfold(args.old_bin, args.input_file, extra_args, args.repeats, args.warmups)
    new_run = run_kinfold(args.new_bin, args.input_file, extra_args, args.repeats, args.warmups)

    summary = {
        "input_file": args.input_file,
        "extra_args": extra_args,
        "repeats": args.repeats,
        "warmups": args.warmups,
        "stdout_equal": old_run["stdout"] == new_run["stdout"],
        "log_equal": old_run["log"] == new_run["log"],
        "stderr_equal": old_run["stderr"] == new_run["stderr"],
        "marker_equal": old_run["parsed"]["markers"] == new_run["parsed"]["markers"],
        "old_median_ns": old_run["median_ns"],
        "new_median_ns": new_run["median_ns"],
        "runtime_ratio_new_over_old": (
            float(new_run["median_ns"]) / float(old_run["median_ns"])
            if old_run["median_ns"]
            else None
        ),
        "old_markers": old_run["parsed"]["markers"],
        "new_markers": new_run["parsed"]["markers"],
    }

    pathlib.Path(args.output).write_text(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

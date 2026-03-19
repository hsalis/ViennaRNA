#!/usr/bin/env python3

import argparse
import json
import subprocess
from pathlib import Path


def run(cmd):
    subprocess.run(cmd, check=True)


def load_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--new-module-dir", required=True)
    parser.add_argument("--old-module-dir", required=True)
    parser.add_argument("--new-worker", required=True)
    parser.add_argument("--old-worker", required=True)
    parser.add_argument("--new-rnafold-bin", required=True)
    parser.add_argument("--old-rnafold-bin", required=True)
    parser.add_argument("--new-rnacofold-bin", required=True)
    parser.add_argument("--old-rnacofold-bin", required=True)
    parser.add_argument("--new-rnasubopt-bin", required=True)
    parser.add_argument("--old-rnasubopt-bin", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--warmups", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=10)
    parser.add_argument("--bootstrap-iterations", type=int, default=5000)
    parser.add_argument("--significance-threshold", type=float, default=0.95)
    parser.add_argument("--expect-avx2", action="store_true")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    python_out = out_dir / "python"
    capi_out = out_dir / "capi"
    cli_out = out_dir / "cli"
    python_out.mkdir(exist_ok=True)
    capi_out.mkdir(exist_ok=True)
    cli_out.mkdir(exist_ok=True)

    run([
        "python3",
        str(Path(__file__).with_name("python_viennarna_compare.py")),
        "--new-module-dir", args.new_module_dir,
        "--old-module-dir", args.old_module_dir,
        "--out-dir", str(python_out),
        "--warmups", str(args.warmups),
        "--repeats", str(args.repeats),
        "--bootstrap-iterations", str(args.bootstrap_iterations),
        "--significance-threshold", str(args.significance_threshold),
    ])

    capi_cmd = [
        "python3",
        str(Path(__file__).with_name("vrna_capi_compare.py")),
        "--new-worker", args.new_worker,
        "--old-worker", args.old_worker,
        "--out-dir", str(capi_out),
        "--warmups", str(args.warmups),
        "--repeats", str(args.repeats),
        "--bootstrap-iterations", str(args.bootstrap_iterations),
        "--significance-threshold", str(args.significance_threshold),
    ]
    if args.expect_avx2:
        capi_cmd.append("--expect-avx2")
    run(capi_cmd)

    run([
        "python3",
        str(Path(__file__).with_name("vrna_cli_compare.py")),
        "--new-rnafold-bin", args.new_rnafold_bin,
        "--old-rnafold-bin", args.old_rnafold_bin,
        "--new-rnacofold-bin", args.new_rnacofold_bin,
        "--old-rnacofold-bin", args.old_rnacofold_bin,
        "--new-rnasubopt-bin", args.new_rnasubopt_bin,
        "--old-rnasubopt-bin", args.old_rnasubopt_bin,
        "--out-dir", str(cli_out),
        "--warmups", str(args.warmups),
        "--repeats", str(args.repeats),
        "--bootstrap-iterations", str(args.bootstrap_iterations),
        "--significance-threshold", str(args.significance_threshold),
    ])

    summary = {
        "python": load_json(python_out / "python_viennarna_compare.summary.json"),
        "capi": load_json(capi_out / "capi_compare.summary.json"),
        "cli": load_json(cli_out / "cli_compare.summary.json"),
    }

    (out_dir / "exactness_perf_gate.summary.json").write_text(json.dumps(summary, indent=2),
                                                              encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

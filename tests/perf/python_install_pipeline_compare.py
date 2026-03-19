#!/usr/bin/env python3

import argparse
import json
import re
import subprocess
from pathlib import Path


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def extract_install_contract(makefile_am: Path):
    text = read_text(makefile_am)
    contract = {}
    for key in (
        "pkgpyexecdir",
        "pkgpycachedir",
        "pkgpyexec_LTLIBRARIES",
        "pkgpyexec_DATA",
        "pkgpycache_DATA",
        "pkgpyvrnaexecdir",
        "pkgpyvrnaexec_DATA",
    ):
        match = re.search(rf"^{re.escape(key)}\s*=\s*(.*?)(?=^\S|\Z)", text, re.M | re.S)
        if match:
            value = re.sub(r"\\\n", " ", match.group(1))
            value = re.sub(r"\s+", " ", value).strip()
            contract[key] = value
    contract["has_install_data_hook"] = "install-data-hook:" in text
    contract["removes__RNA_la"] = "rm -f $(DESTDIR)$(pkgpyexecdir)/_RNA.la" in text
    return contract


def normalize_install_output(output: str) -> list[str]:
    lines = []
    for line in output.splitlines():
        line = line.strip()
        if not line:
            continue
        line = re.sub(r"/opt/venv/lib/python[0-9.]+/site-packages", "$SITE", line)
        line = re.sub(r"/usr/local/lib/python[0-9.]+/site-packages", "$SITE", line)
        line = re.sub(r"/Library/Frameworks/Python.framework/Versions/[0-9.]+/lib/python[0-9.]+/site-packages", "$SITE", line)
        line = re.sub(r"\.\./\.\./\.\./ViennaRNA/interfaces/Python/", "$SRC/", line)
        line = re.sub(r"^make(\[[0-9]+\])?: Entering directory .*$", "", line)
        line = re.sub(r"^make(\[[0-9]+\])?: Leaving directory .*$", "", line)
        if line:
            lines.append(line)
    return lines


def make_install_dry_run(directory: Path):
    proc = subprocess.run(
        ["make", "-n", "-C", str(directory), "install-data"],
        check=True,
        capture_output=True,
        text=True,
    )
    return normalize_install_output(proc.stdout)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--new-dir", default="/Users/howardsalis/ViennaRNA/interfaces/Python")
    parser.add_argument("--old-dir", default="/Users/howardsalis/ViennaRNA/old/ViennaRNA/interfaces/Python")
    parser.add_argument("--out", default="/Users/howardsalis/ViennaRNA/tests/perf/python_install_pipeline_compare.json")
    args = parser.parse_args()

    new_dir = Path(args.new_dir)
    old_dir = Path(args.old_dir)

    new_contract = extract_install_contract(new_dir / "Makefile.am")
    old_contract = extract_install_contract(old_dir / "Makefile.am")

    result = {
        "new_contract": new_contract,
        "old_contract": old_contract,
        "contract_match": new_contract == old_contract,
    }

    try:
        new_install = make_install_dry_run(new_dir)
        result["new_install_dry_run"] = new_install
    except Exception as exc:
        result["new_install_dry_run_error"] = str(exc)
        new_install = None

    try:
        old_install = make_install_dry_run(old_dir)
        result["old_install_dry_run"] = old_install
    except Exception as exc:
        result["old_install_dry_run_error"] = str(exc)
        old_install = None

    if (new_install is not None) and (old_install is not None):
        result["dry_run_match"] = new_install == old_install

    out_path = Path(args.out)
    out_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

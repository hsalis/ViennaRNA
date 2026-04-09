#!/usr/bin/env python3

import argparse
import json
import math
import pathlib
import subprocess
import tempfile


TXSTATE_PREFIX = "#TXSTATE "
TXNEIGH_PREFIX = "#TXNEIGH "
GASCONST = 1.98717
K0 = 273.15


CASES = [
    {
        "name": "elongation_only",
        "sequence": "AUGCUA",
        "args": [
            "--transcription_elongation_rate=1000000",
            "--kinfold_seconds_per_time_unit=1e-5",
            "--max_bubble_width=0",
            "--time=0.5",
            "--num=1",
            "--verbose",
            "--seed",
            "1=2=3",
            "--dump_transcription_neighbors",
        ],
        "expected_neighbor_kinds": {"t"},
        "expected_selected_moves": set(),
        "max_bubble_width": 0,
        "temperature_c": 37.0,
        "phi": 1.0,
        "metropolis": False,
        "elongation_rate_nt_s": 1000000.0,
        "seconds_per_time_unit": 1e-5,
    },
    {
        "name": "hybrid_toggle",
        "sequence": "GGGGAAAAACCCC",
        "args": [
            "--transcription_elongation_rate=100000",
            "--kinfold_seconds_per_time_unit=1e-5",
            "--max_bubble_width=3",
            "--time=10",
            "--num=1",
            "--verbose",
            "--seed",
            "1=2=3",
            "--dump_transcription_neighbors",
        ],
        "expected_neighbor_kinds": {"t", "h", "H"},
        "expected_selected_moves": {"t", "h", "H"},
        "max_bubble_width": 3,
        "temperature_c": 37.0,
        "phi": 1.0,
        "metropolis": False,
        "elongation_rate_nt_s": 100000.0,
        "seconds_per_time_unit": 1e-5,
    },
    {
        "name": "bubble_growth_invasion",
        "sequence": "GGGGAAAAACCCC",
        "args": [
            "--transcription_elongation_rate=100000",
            "--kinfold_seconds_per_time_unit=1e-5",
            "--max_bubble_width=2",
            "--time=30",
            "--num=1",
            "--verbose",
            "--seed",
            "1=2=3",
            "--dump_transcription_neighbors",
        ],
        "expected_neighbor_kinds": {"t", "h", "H", "b", "B"},
        "expected_selected_moves": {"t", "h", "H", "b", "B"},
        "max_bubble_width": 2,
        "temperature_c": 37.0,
        "phi": 1.0,
        "metropolis": False,
        "elongation_rate_nt_s": 100000.0,
        "seconds_per_time_unit": 1e-5,
    },
]


def parse_kv_fields(text):
    result = {}
    for token in text.split():
        key, value = token.split("=", 1)
        result[key] = value
    return result


def parse_trace(stdout, stderr):
    states = []
    current = None
    selected_moves = []

    for line in stderr.splitlines():
        if line.startswith(TXSTATE_PREFIX):
            fields = parse_kv_fields(line[len(TXSTATE_PREFIX):])
            current = {
                "step": int(fields["step"]),
                "tx_len": int(fields["tx_len"]),
                "full_len": int(fields["full_len"]),
                "bubble_left": int(fields["bubble_left"]),
                "bubble_width": int(fields["bubble_width"]),
                "hybrid_left": int(fields["hybrid_left"]),
                "hybrid_len": int(fields["hybrid_len"]),
                "curr_rna_E": float(fields["curr_rna_E"]),
                "curr_total_E": float(fields["curr_total_E"]),
                "seq": fields["seq"],
                "struct": fields["struct"],
                "neighbors": [],
            }
            states.append(current)
        elif line.startswith(TXNEIGH_PREFIX):
            if current is None:
                raise AssertionError("Encountered TXNEIGH before TXSTATE")
            fields = parse_kv_fields(line[len(TXNEIGH_PREFIX):])
            current["neighbors"].append(
                {
                    "kind": fields["kind"],
                    "i": int(fields["i"]),
                    "j": int(fields["j"]),
                    "next_tx_len": int(fields["next_tx_len"]),
                    "next_bubble_left": int(fields["next_bubble_left"]),
                    "next_bubble_width": int(fields["next_bubble_width"]),
                    "next_hybrid_left": int(fields["next_hybrid_left"]),
                    "next_hybrid_len": int(fields["next_hybrid_len"]),
                    "next_rna_E": float(fields["next_rna_E"]),
                    "next_total_E": float(fields["next_total_E"]),
                    "rate": float(fields["rate"]),
                }
            )

    for line in stdout.splitlines():
        parts = line.split()
        if len(parts) >= 6:
            selected_moves.append(parts[4])

    return {
        "states": states,
        "selected_moves": selected_moves,
    }


def run_case(binary, case):
    with tempfile.TemporaryDirectory(prefix=f"kinfold-tx-{case['name']}-") as tmpdir:
        input_path = pathlib.Path(tmpdir) / "case.in"
        input_path.write_text(case["sequence"] + "\n")
        proc = subprocess.run(
            [binary, "--log", str(pathlib.Path(tmpdir) / "kinout")] + case["args"],
            stdin=input_path.open("rb"),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True,
            cwd=tmpdir,
        )
        stdout = proc.stdout.decode("utf-8")
        stderr = proc.stderr.decode("utf-8")
        return {
            "stdout": stdout,
            "stderr": stderr,
            "trace": parse_trace(stdout, stderr),
        }


def transition_rate(case, curr_total, next_total, kind):
    if kind == "t":
        return case["elongation_rate_nt_s"] * case["seconds_per_time_unit"]
    rt = ((case["temperature_c"] + K0) * GASCONST) / 1000.0
    d_e = next_total - curr_total
    if case["metropolis"]:
        return 1.0 if d_e < 0.0 else math.exp(-(d_e / rt * case["phi"]))
    return math.exp(-0.5 * (d_e / rt * case["phi"]))


def validate_bounds(state, neigh):
    if not (0 <= neigh["next_tx_len"] <= state["full_len"]):
        raise AssertionError(f"Invalid next_tx_len for {neigh}")
    if neigh["next_bubble_width"] < 0 or neigh["next_hybrid_len"] < 0:
        raise AssertionError(f"Negative bubble/hybrid size for {neigh}")
    if neigh["next_bubble_width"] > 0:
        if neigh["next_bubble_left"] < 0:
            raise AssertionError(f"Negative bubble_left for {neigh}")
        if neigh["next_bubble_left"] + neigh["next_bubble_width"] > neigh["next_tx_len"]:
            raise AssertionError(f"Bubble exceeds transcript for {neigh}")
    if neigh["next_hybrid_len"] > 0:
        if neigh["next_hybrid_left"] < 0:
            raise AssertionError(f"Negative hybrid_left for {neigh}")
        if neigh["next_hybrid_left"] + neigh["next_hybrid_len"] > neigh["next_tx_len"]:
            raise AssertionError(f"Hybrid exceeds transcript for {neigh}")
        if neigh["next_bubble_width"] > 0:
            if neigh["next_hybrid_left"] < neigh["next_bubble_left"]:
                raise AssertionError(f"Hybrid escapes bubble on left for {neigh}")
            if (neigh["next_hybrid_left"] + neigh["next_hybrid_len"]) > (
                neigh["next_bubble_left"] + neigh["next_bubble_width"]
            ):
                raise AssertionError(f"Hybrid escapes bubble on right for {neigh}")


def validate_neighbor(case, state, neigh):
    validate_bounds(state, neigh)

    kind = neigh["kind"]
    if kind == "t":
        if neigh["next_tx_len"] != state["tx_len"] + 1:
            raise AssertionError(f"Elongation did not advance transcript: {state} -> {neigh}")
        if case["max_bubble_width"] == 0:
            if neigh["next_bubble_width"] != 0 or neigh["next_hybrid_len"] != 0:
                raise AssertionError(f"Non-bubble elongation should not create hybrid state: {neigh}")
        else:
            if neigh["next_hybrid_len"] <= 0:
                raise AssertionError(f"Bubble-mode elongation must pair the new nt: {neigh}")
            if neigh["next_hybrid_left"] + neigh["next_hybrid_len"] != neigh["next_tx_len"]:
                raise AssertionError(f"Hybrid must end at transcript tip after elongation: {neigh}")
            if neigh["next_bubble_width"] > case["max_bubble_width"]:
                raise AssertionError(f"Bubble exceeded width cap: {neigh}")
    elif kind == "h":
        if not (state["bubble_left"] < state["hybrid_left"]):
            raise AssertionError(f"Hybrid form requires exposed upstream bubble nt: {state}")
        if neigh["next_tx_len"] != state["tx_len"]:
            raise AssertionError(f"Hybrid form changed transcript length: {neigh}")
        if neigh["next_hybrid_left"] != state["hybrid_left"] - 1:
            raise AssertionError(f"Hybrid form did not extend upstream by one nt: {neigh}")
        if neigh["next_hybrid_len"] != state["hybrid_len"] + 1:
            raise AssertionError(f"Hybrid form did not grow hybrid length: {neigh}")
    elif kind == "H":
        if state["hybrid_len"] <= 0:
            raise AssertionError(f"Hybrid break requires existing hybrid: {state}")
        if neigh["next_tx_len"] != state["tx_len"]:
            raise AssertionError(f"Hybrid break changed transcript length: {neigh}")
        if neigh["next_hybrid_len"] == 0:
            pass
        else:
            if neigh["next_hybrid_left"] != state["hybrid_left"] + 1:
                raise AssertionError(f"Hybrid break did not retreat upstream edge: {neigh}")
            if neigh["next_hybrid_len"] != state["hybrid_len"] - 1:
                raise AssertionError(f"Hybrid break did not shrink hybrid length: {neigh}")
    elif kind == "b":
        if state["bubble_width"] >= case["max_bubble_width"]:
            raise AssertionError(f"Bubble growth offered at width cap: {state}")
        if state["bubble_left"] <= 0:
            raise AssertionError(f"Bubble growth requires upstream nt: {state}")
        if neigh["next_bubble_left"] != state["bubble_left"] - 1:
            raise AssertionError(f"Bubble growth did not move upstream edge left: {neigh}")
        if neigh["next_bubble_width"] != state["bubble_width"] + 1:
            raise AssertionError(f"Bubble growth did not widen bubble: {neigh}")
        if neigh["next_hybrid_left"] != state["hybrid_left"] - 1:
            raise AssertionError(f"Bubble growth did not extend hybrid upstream: {neigh}")
        if neigh["next_hybrid_len"] != state["hybrid_len"] + 1:
            raise AssertionError(f"Bubble growth did not widen hybrid: {neigh}")
    elif kind == "B":
        if not (state["bubble_left"] < state["hybrid_left"]):
            raise AssertionError(f"Bubble invasion requires exposed bubble nt: {state}")
        if neigh["i"] <= 0 or neigh["j"] <= 0 or neigh["i"] == neigh["j"]:
            raise AssertionError(f"Bubble invasion requires a distinct RNA-RNA pair: {neigh}")
        if neigh["next_tx_len"] != state["tx_len"]:
            raise AssertionError(f"Bubble invasion changed transcript length: {neigh}")
        expected_width = max(state["bubble_width"] - 1, 0)
        if neigh["next_bubble_width"] != expected_width:
            raise AssertionError(f"Bubble invasion did not shrink bubble width correctly: {neigh}")
        if expected_width > 0:
            if neigh["next_bubble_left"] != state["bubble_left"] + 1:
                raise AssertionError(f"Bubble invasion did not advance bubble left edge: {neigh}")
    else:
        raise AssertionError(f"Unexpected special move kind: {kind}")

    expected_rate = transition_rate(case, state["curr_total_E"], neigh["next_total_E"], kind)
    if not math.isclose(neigh["rate"], expected_rate, rel_tol=1e-7, abs_tol=1e-9):
        raise AssertionError(
            f"Rate mismatch for {kind}: observed {neigh['rate']}, expected {expected_rate}"
        )


def validate_case(case, run):
    trace = run["trace"]
    states = trace["states"]
    if not states:
        raise AssertionError(f"No transcription debug states found for case {case['name']}")

    observed_neighbor_kinds = set()
    for state in states:
        for neigh in state["neighbors"]:
            observed_neighbor_kinds.add(neigh["kind"])
            validate_neighbor(case, state, neigh)

    if not case["expected_neighbor_kinds"].issubset(observed_neighbor_kinds):
        raise AssertionError(
            f"Missing expected neighbor kinds for {case['name']}: "
            f"expected {sorted(case['expected_neighbor_kinds'])}, observed {sorted(observed_neighbor_kinds)}"
        )

    selected = {move for move in trace["selected_moves"] if move in {"t", "h", "H", "b", "B"}}
    if not case["expected_selected_moves"].issubset(selected):
        raise AssertionError(
            f"Missing expected selected move kinds for {case['name']}: "
            f"expected {sorted(case['expected_selected_moves'])}, observed {sorted(selected)}"
        )

    return {
        "states": len(states),
        "observed_neighbor_kinds": sorted(observed_neighbor_kinds),
        "observed_selected_moves": sorted(selected),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    summary = {
        "binary": args.binary,
        "cases": [],
        "all_passed": True,
    }

    for case in CASES:
        first = run_case(args.binary, case)
        second = run_case(args.binary, case)
        if first["stdout"] != second["stdout"] or first["stderr"] != second["stderr"]:
            raise AssertionError(f"Seeded determinism failed for case {case['name']}")

        case_summary = validate_case(case, first)
        case_summary["name"] = case["name"]
        case_summary["stdout_lines"] = len([line for line in first["stdout"].splitlines() if line.strip()])
        case_summary["stderr_lines"] = len([line for line in first["stderr"].splitlines() if line.strip()])
        summary["cases"].append(case_summary)

    pathlib.Path(args.output).write_text(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

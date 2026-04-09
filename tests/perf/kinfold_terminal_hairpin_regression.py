#!/usr/bin/env python3

import argparse
import json
import math
import pathlib
import subprocess
import tempfile


HYBRID_ENERGIES = {
    "TT": -0.7,
    "TG": -1.2,
    "TC": -1.5,
    "TA": -0.5,
    "GT": -1.5,
    "GG": -1.7,
    "GC": -2.0,
    "GA": -1.4,
    "CT": -1.3,
    "CG": -1.4,
    "CC": -2.3,
    "CA": -1.6,
    "AT": -0.4,
    "AG": -0.4,
    "AC": -1.4,
    "AA": 0.2,
}

PAIRABLE = {
    ("A", "U"),
    ("U", "A"),
    ("G", "C"),
    ("C", "G"),
    ("G", "U"),
    ("U", "G"),
}


def round_dcal(value):
    return math.copysign(math.floor(abs(value) * 100.0 + 0.5), value) / 100.0


def build_pair_tables(structure):
    pair_table = [-1] * len(structure)
    parent_open = [-1] * len(structure)
    stack = []
    for idx, char in enumerate(structure):
        if char == "(":
            parent_open[idx] = stack[-1] if stack else -1
            stack.append(idx)
        elif char == ")":
            if not stack:
                raise AssertionError(f"Unbalanced structure: {structure}")
            left = stack.pop()
            pair_table[left] = idx
            pair_table[idx] = left
            parent_open[idx] = parent_open[left]
    if stack:
        raise AssertionError(f"Unbalanced structure: {structure}")
    return pair_table, parent_open


def find_terminal_component(structure, bubble_left):
    pair_table, parent_open = build_pair_tables(structure)
    best = None
    for idx, mate in enumerate(pair_table):
        if mate <= idx:
            continue
        if mate >= bubble_left:
            continue
        if parent_open[idx] != -1:
            continue
        if best is None or mate > best[1]:
            best = (idx, mate)
    return best


def eval_substructure(RNA, sequence, structure):
    fc = RNA.fold_compound(sequence.replace("T", "U"))
    return round_dcal(fc.eval_structure(structure))


def terminal_hairpin_energy(RNA, sequence, structure, bubble_left):
    root = find_terminal_component(structure, bubble_left)
    if root is None:
        return 0.0, None
    left, right = root
    return (
        eval_substructure(RNA, sequence[left:right + 1], structure[left:right + 1]),
        root,
    )


def decorate_structure(structure, bubble_left, bubble_length):
    chars = list(structure)
    if bubble_length > 0:
        chars[bubble_left] = "b"
        chars[bubble_left + bubble_length - 1] = "b"
    return "".join(chars)


def template_dna(full_sequence):
    mapping = {"A": "T", "U": "A", "G": "C", "C": "G", "T": "A"}
    return "".join(mapping[nt] for nt in full_sequence.upper())


def hybrid_energy(full_sequence, hybrid_left, hybrid_len):
    if hybrid_len <= 1:
        return 0.0
    template = template_dna(full_sequence)
    total = 0.0
    for idx in range(hybrid_left, hybrid_left + hybrid_len - 1):
        total += HYBRID_ENERGIES[template[idx:idx + 2]]
    return round_dcal(total)


def smoothed_derivative(energies, times, time_step):
    tau = 5.0 * time_step
    values = []
    smoothed = 0.0
    previous = 0.0
    initialized = False
    previous_time = 0.0
    for energy, sample_time in zip(energies, times):
        if not initialized:
            smoothed = energy
            previous = smoothed
            values.append(0.0)
            initialized = True
        else:
            dt = sample_time - previous_time
            if dt <= 0.0:
                values.append(0.0)
            else:
                alpha = 1.0 - math.exp(-dt / tau)
                smoothed = alpha * energy + (1.0 - alpha) * smoothed
                values.append((smoothed - previous) / dt)
                previous = smoothed
        previous_time = sample_time
    return values


def is_pairable(left_nt, right_nt):
    return (left_nt.upper().replace("T", "U"), right_nt.upper().replace("T", "U")) in PAIRABLE


def plus1_is_valid(sequence, structure, bubble_left, hybrid_left, hybrid_len, root):
    if root is None or hybrid_len <= 0 or hybrid_left < 0:
        return False
    left, _right = root
    if left <= 0 or (left - 1) >= bubble_left:
        return False
    pair_table, _parent_open = build_pair_tables(structure)
    if pair_table[left - 1] != -1:
        return False
    if structure[hybrid_left] != ".":
        return False
    if abs(hybrid_left - (left - 1)) < 1:
        return False
    return is_pairable(sequence[left - 1], sequence[hybrid_left])


def parse_cli_line(line):
    parts = line.split()
    if len(parts) < 10:
        raise AssertionError(f"Unexpected CLI line: {line}")
    return {
        "structure": parts[0],
        "energy": float(parts[1]),
        "rna_dna": float(parts[2]),
        "dna_dna": float(parts[3]),
        "terminal_hairpin": float(parts[4]),
        "terminal_hairpin_dt": float(parts[5]),
        "terminal_hairpin_plus1": float(parts[6]),
        "rna_dna_minus1_ddg": float(parts[7]),
        "plus1_ddg": float(parts[8]),
        "time": float(parts[9]),
        "marker": parts[10] if len(parts) > 10 else "",
    }


def run_cli(binary, sequence, extra_args):
    with tempfile.TemporaryDirectory(prefix="kinfold-terminal-hairpin-") as tmpdir:
        proc = subprocess.run(
            [binary] + list(extra_args),
            input=(sequence + "\n").encode("utf-8"),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            cwd=tmpdir,
            check=True,
        )
        stdout = proc.stdout.decode("utf-8")
        stderr = proc.stderr.decode("utf-8")
        lines = [parse_cli_line(line) for line in stdout.splitlines() if line.strip()]
        return {"stdout": stdout, "stderr": stderr, "lines": lines}


def assert_close(observed, expected, label, abs_tol=1e-6):
    if not math.isclose(observed, expected, rel_tol=1e-6, abs_tol=abs_tol):
        raise AssertionError(f"{label}: observed {observed}, expected {expected}")


def validate_python_metrics(RNA, traj):
    sequence = traj["sequence"]
    time_step = traj["parameters"]["time_step"]

    expected_terminal = []
    roots = []
    for structure, bubble_left in zip(traj["structure"], traj["bubble_left"]):
        value, root = terminal_hairpin_energy(RNA, sequence, structure.replace("x", "."), bubble_left)
        expected_terminal.append(value)
        roots.append(root)

    expected_dt = smoothed_derivative(expected_terminal, traj["time"], time_step)
    expected_minus1 = []
    for idx in range(len(traj["time"])):
        current_hybrid = round_dcal(traj["total_energy"][idx] - traj["energy"][idx])
        hybrid_len = traj["hybrid_length"][idx]
        hybrid_left = traj["hybrid_left"][idx]
        if hybrid_len > 0:
            minus1 = hybrid_energy(sequence, hybrid_left + 1, hybrid_len - 1)
            expected_minus1.append(round_dcal(minus1 - current_hybrid))
        else:
            expected_minus1.append(0.0)

    for idx, (observed, expected) in enumerate(zip(traj["terminal_hairpin_energy"], expected_terminal)):
        assert_close(observed, expected, f"terminal_hairpin_energy[{idx}]", abs_tol=1e-2)

    for idx, (observed, expected) in enumerate(zip(traj["terminal_hairpin_dg_dt"], expected_dt)):
        assert_close(observed, expected, f"terminal_hairpin_dg_dt[{idx}]", abs_tol=1e-6)

    for idx, (observed, expected) in enumerate(zip(traj["rna_dna_minus1_ddg"], expected_minus1)):
        assert_close(observed, expected, f"rna_dna_minus1_ddg[{idx}]", abs_tol=1e-2)

    valid_plus1_samples = 0
    for idx, observed in enumerate(traj["plus1_ddg"]):
        if plus1_is_valid(
            sequence,
            traj["structure"][idx].replace("x", "."),
            traj["bubble_left"][idx],
            traj["hybrid_left"][idx],
            traj["hybrid_length"][idx],
            roots[idx],
        ):
            valid_plus1_samples += 1
            assert_close(
                observed,
                traj["terminal_hairpin_plus1_energy"][idx]
                - traj["terminal_hairpin_energy"][idx]
                - traj["rna_dna_minus1_ddg"][idx],
                f"plus1_ddg[{idx}]",
                abs_tol=1e-2,
            )
        else:
            assert_close(observed, 0.0, f"plus1_ddg[{idx}] invalid", abs_tol=1e-9)

    return {
        "terminal_hairpin_energy_matches": True,
        "terminal_hairpin_dg_dt_matches": True,
        "rna_dna_minus1_ddg_matches": True,
        "plus1_identity_matches": True,
        "nonzero_terminal_samples": sum(1 for value in traj["terminal_hairpin_energy"] if abs(value) > 1e-9),
        "valid_plus1_samples": valid_plus1_samples,
        "nonzero_plus1_samples": sum(1 for value in traj["plus1_ddg"] if abs(value) > 1e-9),
        "roots": roots,
    }


def validate_cli_python_parity(cli_lines, traj):
    if len(cli_lines) != len(traj["time"]) - 1:
        raise AssertionError(
            f"CLI/Python sample count mismatch: CLI {len(cli_lines)} vs Python {len(traj['time']) - 1}"
        )

    for idx, line in enumerate(cli_lines, start=1):
        expected_structure = decorate_structure(
            traj["structure"][idx], traj["bubble_left"][idx], traj["bubble_length"][idx]
        )
        if line["structure"] != expected_structure:
            raise AssertionError(
                f"CLI structure mismatch at sample {idx}: {line['structure']} vs {expected_structure}"
            )
        assert_close(line["time"], traj["time"][idx], f"CLI time[{idx}]", abs_tol=1e-3)
        assert_close(line["energy"], traj["energy"][idx], f"CLI RNA energy[{idx}]", abs_tol=1e-2)
        assert_close(
            line["rna_dna"],
            traj["total_energy"][idx] - traj["energy"][idx],
            f"CLI RNA:DNA energy[{idx}]",
            abs_tol=1e-2,
        )
        assert_close(line["dna_dna"], traj["dna_duplex_energy"][idx], f"CLI DNA:DNA energy[{idx}]", abs_tol=1e-2)
        assert_close(
            line["terminal_hairpin"],
            traj["terminal_hairpin_energy"][idx],
            f"CLI terminal hairpin[{idx}]",
            abs_tol=1e-2,
        )
        assert_close(
            line["terminal_hairpin_dt"],
            traj["terminal_hairpin_dg_dt"][idx],
            f"CLI terminal hairpin dt[{idx}]",
            abs_tol=1e-4,
        )
        assert_close(
            line["terminal_hairpin_plus1"],
            traj["terminal_hairpin_plus1_energy"][idx],
            f"CLI terminal hairpin plus1[{idx}]",
            abs_tol=1e-2,
        )
        assert_close(
            line["rna_dna_minus1_ddg"],
            traj["rna_dna_minus1_ddg"][idx],
            f"CLI RNA:DNA minus1 ddG[{idx}]",
            abs_tol=1e-2,
        )
        assert_close(line["plus1_ddg"], traj["plus1_ddg"][idx], f"CLI plus1 ddG[{idx}]", abs_tol=1e-2)

    return {"cli_python_parity": True, "samples_compared": len(cli_lines)}


def validate_double_hairpin_choice(RNA):
    sequence = "GGGGAAAACCCCGGGGAAAACCCCUUUU"
    traj = RNA.kinfold(
        sequence,
        time=4000.0,
        time_step=400.0,
        seed=(1, 2, 3),
        silent=True,
        transcription_elongation_rate=1000.0,
        max_bubble_width=3,
    )
    final_idx = len(traj["time"]) - 1
    structure = traj["structure"][final_idx].replace("x", ".")
    root = find_terminal_component(structure, traj["bubble_left"][final_idx])
    if root is None:
        raise AssertionError("Expected a downstream terminal hairpin in double_hairpin case")
    left, right = root
    if (left, right) != (12, 23):
        raise AssertionError(f"Unexpected downstream hairpin root: {(left, right)}")
    expected = eval_substructure(RNA, sequence[left:right + 1], structure[left:right + 1])
    assert_close(
        traj["terminal_hairpin_energy"][final_idx],
        expected,
        "double_hairpin terminal energy",
        abs_tol=1e-2,
    )
    return {
        "final_root": [left, right],
        "final_terminal_hairpin_energy": traj["terminal_hairpin_energy"][final_idx],
    }


def validate_nonbubble_zeroes(RNA):
    traj = RNA.kinfold("GGGAAACCCUUU", time=200.0, seed=(1, 2, 3), silent=True)
    for field in (
        "terminal_hairpin_energy",
        "terminal_hairpin_dg_dt",
        "terminal_hairpin_plus1_energy",
        "rna_dna_minus1_ddg",
        "plus1_ddg",
    ):
        if any(abs(value) > 1e-9 for value in traj[field]):
            raise AssertionError(f"Expected zeros for non-bubble field {field}: {traj[field]}")
    return {"nonbubble_zero_fields": True, "samples": len(traj["time"])}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    import RNA

    sequence = "GGGGGGAAAAAACCCCCCUGUGUG"
    shared_args = {
        "time": 2000.0,
        "time_step": 200.0,
        "seed": (1, 2, 3),
        "silent": True,
        "transcription_elongation_rate": 1000.0,
        "max_bubble_width": 3,
    }

    python_traj = RNA.kinfold(sequence, **shared_args)
    cli = run_cli(
        args.binary,
        sequence,
        [
            "--num",
            "1",
            "--time",
            "2000",
            "--time-step",
            "200",
            "--seed",
            "1=2=3",
            "--transcription_elongation_rate",
            "1000",
            "--max_bubble_width",
            "3",
        ],
    )

    summary = {
        "binary": args.binary,
        "python_metric_validation": validate_python_metrics(RNA, python_traj),
        "cli_python_parity": validate_cli_python_parity(cli["lines"], python_traj),
        "double_hairpin_choice": validate_double_hairpin_choice(RNA),
        "nonbubble_zeroes": validate_nonbubble_zeroes(RNA),
        "designed_terminal_last_sample": {
            "time": python_traj["time"][-1],
            "structure": python_traj["structure"][-1],
            "terminal_hairpin_energy": python_traj["terminal_hairpin_energy"][-1],
            "terminal_hairpin_dg_dt": python_traj["terminal_hairpin_dg_dt"][-1],
            "terminal_hairpin_plus1_energy": python_traj["terminal_hairpin_plus1_energy"][-1],
            "rna_dna_minus1_ddg": python_traj["rna_dna_minus1_ddg"][-1],
            "plus1_ddg": python_traj["plus1_ddg"][-1],
        },
    }

    pathlib.Path(args.output).write_text(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

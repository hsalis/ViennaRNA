#!/usr/bin/env python3

import argparse
import json
import os


def main():
    parser = argparse.ArgumentParser(description="Smoke-test the Python Kinfold interface")
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    import RNA

    single = RNA.kinfold(
        "GGGAAACCC",
        time=50.0,
        seed=(1, 2, 3),
        silent=True,
    )

    transcription = RNA.kinfold(
        "GGGAAACCCUUU",
        time=1000.0,
        seed=(1, 2, 3),
        silent=True,
        transcription_elongation_rate=1000.0,
        max_bubble_width=3,
    )

    designed_terminal = RNA.kinfold(
        "GGGGGGAAAAAACCCCCCUGUGUG",
        time=2000.0,
        time_step=200.0,
        seed=(1, 2, 3),
        silent=True,
        transcription_elongation_rate=1000.0,
        max_bubble_width=3,
    )

    batch = RNA.kinfold_batch(
        ["GGGAAACCC", "GGGAAACCCUUU"],
        trajectories_per_sequence=2,
        time=1.0,
        silent=True,
        seeds=[(1, 2, 3), (4, 5, 6)],
    )

    summary = {
        "single": {
            "keys": sorted(single.keys()),
            "termination": single["termination"],
            "step_count": len(single["time"]),
            "moves_head": single["move"][:5],
        },
        "transcription": {
            "termination": transcription["termination"],
            "step_count": len(transcription["time"]),
            "sequence_state_head": transcription["sequence_state"][:5],
            "transcribed_length_head": transcription["transcribed_length"][:5],
            "bubble_length_head": transcription["bubble_length"][:5],
            "terminal_hairpin_energy_head": transcription["terminal_hairpin_energy"][:5],
            "terminal_hairpin_dg_dt_head": transcription["terminal_hairpin_dg_dt"][:5],
            "terminal_hairpin_plus1_energy_head": transcription["terminal_hairpin_plus1_energy"][:5],
            "rna_dna_minus1_ddg_head": transcription["rna_dna_minus1_ddg"][:5],
            "plus1_ddg_head": transcription["plus1_ddg"][:5],
            "move_head": transcription["move"][:10],
        },
        "designed_terminal": {
            "termination": designed_terminal["termination"],
            "step_count": len(designed_terminal["time"]),
            "terminal_hairpin_energy_tail": designed_terminal["terminal_hairpin_energy"][-5:],
            "terminal_hairpin_dg_dt_tail": designed_terminal["terminal_hairpin_dg_dt"][-5:],
            "terminal_hairpin_plus1_energy_tail": designed_terminal["terminal_hairpin_plus1_energy"][-5:],
            "rna_dna_minus1_ddg_tail": designed_terminal["rna_dna_minus1_ddg"][-5:],
            "plus1_ddg_tail": designed_terminal["plus1_ddg"][-5:],
        },
        "batch": {
            "group_count": len(batch),
            "trajectory_counts": [len(group["trajectories"]) for group in batch],
            "terminations": [[traj["termination"] for traj in group["trajectories"]] for group in batch],
        },
        "module_path": getattr(RNA, "__file__", None),
        "cwd": os.getcwd(),
    }

    with open(args.output, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()

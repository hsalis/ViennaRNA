#!/usr/bin/env python3

import argparse
import csv
import json
import pathlib
import random
import statistics
import time


U_TRACT = "UUUUUUUUU"
GC_TRACT = "GCGCGCGCG"
LOOP = "AAAAA"


def make_seed_triplet(rng):
    return tuple(rng.randint(1, 65535) for _ in range(3))


def stem_sequence(length):
    return ("GCGC" * 8)[:length]


def build_cases():
    cases = []
    for stem_len in range(4, 9):
        stem = stem_sequence(stem_len)
        terminator = f"{stem}{LOOP}{stem}{U_TRACT}"
        unstructured_control = f"{'C' * (2 * stem_len + len(LOOP))}{U_TRACT}"
        GC_tract_control = f"{stem}{LOOP}{stem}{GC_TRACT}"
        unstructured_GC_tract_control = f"{'C' * (2 * stem_len + len(LOOP))}{GC_TRACT}"

        cases.append(
            {
                "case_id": f"terminator_stem_{stem_len}",
                "case_type": "terminator",
                "stem_length": stem_len,
                "sequence": terminator,
                "description": f"GC-rich hairpin [length {stem_len}] upstream of canonical U-tract",
            }
        )
        cases.append(
            {
                "case_id": f"control_unstructured_{stem_len}",
                "case_type": "control",
                "stem_length": stem_len,
                "sequence": unstructured_control,
                "description": f"Unstructured RNA region [length {stem_len}] upstream of canonical U-tract",
            }
        )
        cases.append(
            {
                "case_id": f"control_unstructured_{stem_len}",
                "case_type": "control",
                "stem_length": stem_len,
                "sequence": GC_tract_control,
                "description": f"GC-rich hairpin [length {stem_len}] upstream of a GC-rich tract",
            }
        )
        cases.append(
            {
                "case_id": f"control_unstructured_{stem_len}",
                "case_type": "control",
                "stem_length": stem_len,
                "sequence": unstructured_GC_tract_control,
                "description": f"Unstructured RNA region [length {stem_len}] upstream of a GC-rich tract",
            }
        )
    return cases


def safe_mean(values):
    return statistics.fmean(values) if values else None


def trajectory_summary(traj):
    hairpin_times = [
        sample_time
        for sample_time, energy in zip(traj["time"], traj["terminal_hairpin_energy"])
        if energy < 0.0
    ]
    terminal_present_times = [
        sample_time
        for sample_time, energy in zip(traj["time"], traj["terminal_hairpin_energy"])
        if abs(energy) > 1e-9
    ]
    dg_dt = traj["terminal_hairpin_dg_dt"]
    min_dg_dt = min(dg_dt) if dg_dt else 0.0
    min_dg_dt_index = dg_dt.index(min_dg_dt) if dg_dt else 0
    negative_dg_dt = [value for value in dg_dt if value < 0.0]
    return {
        "trajectory_index": traj["trajectory_index"],
        "termination": traj["termination"],
        "final_time": traj["time"][-1] if traj["time"] else 0.0,
        "final_rna_energy": traj["energy"][-1] if traj["energy"] else 0.0,
        "final_total_energy": traj["total_energy"][-1] if traj["total_energy"] else 0.0,
        "final_terminal_hairpin_energy": (
            traj["terminal_hairpin_energy"][-1] if traj["terminal_hairpin_energy"] else 0.0
        ),
        "final_plus1_ddg": traj["plus1_ddg"][-1] if traj["plus1_ddg"] else 0.0,
        "min_rna_energy": min(traj["energy"]) if traj["energy"] else 0.0,
        "min_total_energy": min(traj["total_energy"]) if traj["total_energy"] else 0.0,
        "min_terminal_hairpin_energy": (
            min(traj["terminal_hairpin_energy"]) if traj["terminal_hairpin_energy"] else 0.0
        ),
        "min_plus1_ddg": min(traj["plus1_ddg"]) if traj["plus1_ddg"] else 0.0,
        "hairpin_detected": any(abs(value) > 1e-9 for value in traj["terminal_hairpin_energy"]),
        "hairpin_formed": any(value < 0.0 for value in traj["terminal_hairpin_energy"]),
        "first_terminal_hairpin_time": hairpin_times[0] if hairpin_times else None,
        "first_terminal_component_time": terminal_present_times[0] if terminal_present_times else None,
        "min_terminal_hairpin_dg_dt": min_dg_dt,
        "max_negative_terminal_hairpin_dg_dt_magnitude": abs(min(negative_dg_dt)) if negative_dg_dt else 0.0,
        "time_of_min_terminal_hairpin_dg_dt": (
            traj["time"][min_dg_dt_index] if traj["time"] else None
        ),
    }


def summarize_case(case, trajectories, static_mfe):
    rows = [trajectory_summary(traj) for traj in trajectories]
    hairpin_times = [row["first_terminal_hairpin_time"] for row in rows if row["first_terminal_hairpin_time"] is not None]
    component_times = [
        row["first_terminal_component_time"] for row in rows if row["first_terminal_component_time"] is not None
    ]
    return {
        **case,
        "static_mfe_energy": static_mfe["energy"],
        "static_mfe_structure": static_mfe["structure"],
        "trajectory_count": len(rows),
        "hairpin_detected_count": sum(int(row["hairpin_detected"]) for row in rows),
        "hairpin_formed_count": sum(int(row["hairpin_formed"]) for row in rows),
        "hairpin_detected_rate": (
            float(sum(int(row["hairpin_detected"]) for row in rows)) / float(len(rows)) if rows else None
        ),
        "hairpin_formed_rate": (
            float(sum(int(row["hairpin_formed"]) for row in rows)) / float(len(rows)) if rows else None
        ),
        "mean_first_terminal_component_time": safe_mean(component_times),
        "mean_first_terminal_hairpin_time": safe_mean(hairpin_times),
        "mean_min_rna_energy": safe_mean([row["min_rna_energy"] for row in rows]),
        "best_min_rna_energy": min((row["min_rna_energy"] for row in rows), default=0.0),
        "mean_min_terminal_hairpin_energy": safe_mean(
            [row["min_terminal_hairpin_energy"] for row in rows]
        ),
        "best_min_terminal_hairpin_energy": min(
            (row["min_terminal_hairpin_energy"] for row in rows), default=0.0
        ),
        "mean_min_terminal_hairpin_dg_dt": safe_mean(
            [row["min_terminal_hairpin_dg_dt"] for row in rows]
        ),
        "best_min_terminal_hairpin_dg_dt": min(
            (row["min_terminal_hairpin_dg_dt"] for row in rows), default=0.0
        ),
        "mean_max_negative_terminal_hairpin_dg_dt_magnitude": safe_mean(
            [row["max_negative_terminal_hairpin_dg_dt_magnitude"] for row in rows]
        ),
        "mean_time_of_min_terminal_hairpin_dg_dt": safe_mean(
            [
                row["time_of_min_terminal_hairpin_dg_dt"]
                for row in rows
                if row["time_of_min_terminal_hairpin_dg_dt"] is not None
            ]
        ),
        "mean_final_rna_energy": safe_mean([row["final_rna_energy"] for row in rows]),
        "mean_final_terminal_hairpin_energy": safe_mean(
            [row["final_terminal_hairpin_energy"] for row in rows]
        ),
        "mean_min_plus1_ddg": safe_mean([row["min_plus1_ddg"] for row in rows]),
        "best_min_plus1_ddg": min((row["min_plus1_ddg"] for row in rows), default=0.0),
        "trajectory_summaries": rows,
    }


def summarize_overall(case_summaries):
    grouped = {"terminator": [], "control": []}
    for case in case_summaries:
        grouped[case["case_type"]].append(case)

    overall = {}
    for case_type, rows in grouped.items():
        overall[case_type] = {
            "case_count": len(rows),
            "mean_static_mfe_energy": safe_mean([row["static_mfe_energy"] for row in rows]),
            "mean_hairpin_formed_rate": safe_mean([row["hairpin_formed_rate"] for row in rows]),
            "mean_min_rna_energy": safe_mean([row["mean_min_rna_energy"] for row in rows]),
            "mean_min_terminal_hairpin_energy": safe_mean(
                [row["mean_min_terminal_hairpin_energy"] for row in rows]
            ),
            "mean_min_terminal_hairpin_dg_dt": safe_mean(
                [row["mean_min_terminal_hairpin_dg_dt"] for row in rows]
            ),
            "mean_max_negative_terminal_hairpin_dg_dt_magnitude": safe_mean(
                [row["mean_max_negative_terminal_hairpin_dg_dt_magnitude"] for row in rows]
            ),
            "mean_time_of_min_terminal_hairpin_dg_dt": safe_mean(
                [
                    row["mean_time_of_min_terminal_hairpin_dg_dt"]
                    for row in rows
                    if row["mean_time_of_min_terminal_hairpin_dg_dt"] is not None
                ]
            ),
            "mean_first_terminal_hairpin_time": safe_mean(
                [row["mean_first_terminal_hairpin_time"] for row in rows if row["mean_first_terminal_hairpin_time"] is not None]
            ),
        }

    terminator_rows = grouped["terminator"]
    control_rows = grouped["control"]
    terminator_mean_min = safe_mean([row["mean_min_terminal_hairpin_dg_dt"] for row in terminator_rows])
    control_mean_min = safe_mean([row["mean_min_terminal_hairpin_dg_dt"] for row in control_rows])
    overall["derivative_comparison"] = {
        "terminator_mean_min_terminal_hairpin_dg_dt": terminator_mean_min,
        "control_mean_min_terminal_hairpin_dg_dt": control_mean_min,
        "delta_mean_min_terminal_hairpin_dg_dt_terminator_minus_control": (
            terminator_mean_min - control_mean_min
            if (terminator_mean_min is not None and control_mean_min is not None)
            else None
        ),
        "terminators_more_negative_on_average": (
            terminator_mean_min < control_mean_min
            if (terminator_mean_min is not None and control_mean_min is not None)
            else None
        ),
        "terminator_mean_max_negative_terminal_hairpin_dg_dt_magnitude": safe_mean(
            [row["mean_max_negative_terminal_hairpin_dg_dt_magnitude"] for row in terminator_rows]
        ),
        "control_mean_max_negative_terminal_hairpin_dg_dt_magnitude": safe_mean(
            [row["mean_max_negative_terminal_hairpin_dg_dt_magnitude"] for row in control_rows]
        ),
    }

    return overall


def write_csv(path, case_summaries):
    fieldnames = [
        "case_id",
        "case_type",
        "stem_length",
        "sequence",
        "static_mfe_energy",
        "static_mfe_structure",
        "trajectory_count",
        "hairpin_detected_count",
        "hairpin_formed_count",
        "hairpin_detected_rate",
        "hairpin_formed_rate",
        "mean_first_terminal_component_time",
        "mean_first_terminal_hairpin_time",
        "mean_min_rna_energy",
        "best_min_rna_energy",
        "mean_min_terminal_hairpin_energy",
        "best_min_terminal_hairpin_energy",
        "mean_min_terminal_hairpin_dg_dt",
        "best_min_terminal_hairpin_dg_dt",
        "mean_max_negative_terminal_hairpin_dg_dt_magnitude",
        "mean_time_of_min_terminal_hairpin_dg_dt",
        "mean_final_rna_energy",
        "mean_final_terminal_hairpin_energy",
        "mean_min_plus1_ddg",
        "best_min_plus1_ddg",
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for case in case_summaries:
            writer.writerow({name: case.get(name) for name in fieldnames})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--output-csv")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--trajectories-per-sequence", type=int, default=25)
    parser.add_argument("--time", type=float, default=10000.0)
    parser.add_argument("--time-step", type=float, default=1.0)
    parser.add_argument("--transcription-elongation-rate", type=float, default=3000.0)
    parser.add_argument("--kinfold-seconds-per-time-unit", type=float, default=1e-5)
    parser.add_argument("--max-bubble-width", type=int, default=9)
    parser.add_argument("--temperature", type=float, default=37.0)
    parser.add_argument("--dangle", type=int, default=2)
    args = parser.parse_args()

    if args.trajectories_per_sequence < 1:
        raise ValueError("--trajectories-per-sequence must be >= 1")

    import RNA

    cases = build_cases()
    rng = random.Random(args.seed)
    seeds = [make_seed_triplet(rng) for _ in cases]
    sequences = [case["sequence"] for case in cases]

    static_mfe = {}
    for case in cases:
        fc = RNA.fold_compound(case["sequence"])
        structure, energy = fc.mfe()
        static_mfe[case["case_id"]] = {"structure": structure, "energy": energy}

    start_ns = time.perf_counter_ns()
    groups = RNA.kinfold_batch(
        sequences,
        trajectories_per_sequence=args.trajectories_per_sequence,
        time=args.time,
        time_step=args.time_step,
        Temp=args.temperature,
        dangle=args.dangle,
        silent=True,
        transcription_elongation_rate=args.transcription_elongation_rate,
        kinfold_seconds_per_time_unit=args.kinfold_seconds_per_time_unit,
        max_bubble_width=args.max_bubble_width,
        seeds=seeds,
    )
    wall_ns = time.perf_counter_ns() - start_ns

    case_summaries = []
    for case, group in zip(cases, groups):
        case_summaries.append(
            summarize_case(case, group["trajectories"], static_mfe[case["case_id"]])
        )

    payload = {
        "module_path": getattr(RNA, "__file__", None),
        "seed": args.seed,
        "trajectories_per_sequence": args.trajectories_per_sequence,
        "time": args.time,
        "time_step": args.time_step,
        "transcription_elongation_rate": args.transcription_elongation_rate,
        "kinfold_seconds_per_time_unit": args.kinfold_seconds_per_time_unit,
        "max_bubble_width": args.max_bubble_width,
        "temperature": args.temperature,
        "dangle": args.dangle,
        "wall_ns": wall_ns,
        "case_count": len(case_summaries),
        "cases": case_summaries,
        "overall": summarize_overall(case_summaries),
    }

    pathlib.Path(args.output_json).write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    if args.output_csv:
        write_csv(args.output_csv, case_summaries)


if __name__ == "__main__":
    main()

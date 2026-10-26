#!/usr/bin/env python3

import argparse
import csv
import json
import math
import pathlib
import random
import statistics
import subprocess
import time

import numpy as np
from scipy.stats import pearsonr, spearmanr
from sklearn.linear_model import Lasso, LassoCV
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import KFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from kinfold_intrinsic_terminator_benchmark import make_seed_triplet, summarize_case


DEFAULT_FEATURES = [
    "static_mfe_energy",
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


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Run Kinfold transcription-aware intrinsic terminator metrics on a table of sequences "
            "with measured termination efficiencies, then fit a standardized LASSO regression model."
        )
    )
    parser.add_argument("--input-csv", required=True, help="CSV with at least sequence and efficiency columns")
    parser.add_argument("--output-json", required=True)
    parser.add_argument("--output-csv", help="Optional per-sequence feature/prediction CSV")
    parser.add_argument("--id-column", default="id")
    parser.add_argument("--sequence-column", default="sequence")
    parser.add_argument("--efficiency-column", default="efficiency")
    parser.add_argument("--feature-columns", nargs="+", default=DEFAULT_FEATURES)
    parser.add_argument("--backend", choices=("auto", "python", "cli"), default="auto")
    parser.add_argument("--binary", help="Path to Kinfold CLI binary for cli backend or auto fallback")
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--trajectories-per-sequence", type=int, default=50)
    parser.add_argument("--time", type=float, default=10000.0)
    parser.add_argument("--time-step", type=float, default=1.0)
    parser.add_argument("--transcription-elongation-rate", type=float, default=3000.0)
    parser.add_argument("--kinfold-seconds-per-time-unit", type=float, default=1e-5)
    parser.add_argument("--max-bubble-width", type=int, default=9)
    parser.add_argument("--temperature", type=float, default=37.0)
    parser.add_argument("--dangle", type=int, default=2)
    parser.add_argument("--cv-folds", type=int, default=5)
    parser.add_argument("--max-iter", type=int, default=20000)
    return parser.parse_args()


def safe_mean(values):
    return statistics.fmean(values) if values else None


def load_rows(path, id_column, sequence_column, efficiency_column):
    rows = []
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if sequence_column not in reader.fieldnames:
            raise ValueError(f"Missing sequence column '{sequence_column}' in {path}")
        if efficiency_column not in reader.fieldnames:
            raise ValueError(f"Missing efficiency column '{efficiency_column}' in {path}")

        for idx, row in enumerate(reader, start=1):
            sequence = row[sequence_column].strip().upper().replace("T", "U")
            if not sequence:
                raise ValueError(f"Empty sequence in row {idx}")
            try:
                efficiency = float(row[efficiency_column])
            except Exception as exc:
                raise ValueError(f"Non-numeric efficiency in row {idx}: {row[efficiency_column]!r}") from exc

            rows.append(
                {
                    "id": row.get(id_column) or f"row_{idx:03d}",
                    "sequence": sequence,
                    "measured_efficiency": efficiency,
                    "source_row": idx,
                }
            )
    if len(rows) < 2:
        raise ValueError("Need at least 2 sequences for regression")
    return rows


def compute_static_mfe(RNA, sequence):
    fc = RNA.fold_compound(sequence)
    structure, energy = fc.mfe()
    return {"structure": structure, "energy": float(energy)}


def parse_cli_trajectories(stdout):
    trajectories = []
    current = {
        "time": [],
        "energy": [],
        "total_energy": [],
        "terminal_hairpin_energy": [],
        "terminal_hairpin_dg_dt": [],
        "terminal_hairpin_plus1_energy": [],
        "plus1_ddg": [],
        "trajectory_index": 0,
    }

    for raw_line in stdout.splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 10:
            continue

        rna_energy = float(parts[1])
        rna_dna = float(parts[2])
        current["time"].append(float(parts[9]))
        current["energy"].append(rna_energy)
        current["total_energy"].append(rna_energy + rna_dna)
        current["terminal_hairpin_energy"].append(float(parts[4]))
        current["terminal_hairpin_dg_dt"].append(float(parts[5]))
        current["terminal_hairpin_plus1_energy"].append(float(parts[6]))
        current["plus1_ddg"].append(float(parts[8]))
        marker = parts[10] if len(parts) > 10 else ""
        if marker:
            current["termination"] = marker
            trajectories.append(current)
            current = {
                "time": [],
                "energy": [],
                "total_energy": [],
                "terminal_hairpin_energy": [],
                "terminal_hairpin_dg_dt": [],
                "terminal_hairpin_plus1_energy": [],
                "plus1_ddg": [],
                "trajectory_index": len(trajectories),
            }

    return trajectories


def run_kinfold_cli(binary, sequence, seed_triplet, args):
    cmd = [
        binary,
        "--num",
        str(args.trajectories_per_sequence),
        "--time",
        str(args.time),
        "--time-step",
        str(args.time_step),
        "--seed",
        f"{seed_triplet[0]}={seed_triplet[1]}={seed_triplet[2]}",
        "--transcription_elongation_rate",
        str(args.transcription_elongation_rate),
        "--kinfold_seconds_per_time_unit",
        str(args.kinfold_seconds_per_time_unit),
        "--max_bubble_width",
        str(args.max_bubble_width),
        "--Temp",
        str(args.temperature),
        "--dangle",
        str(args.dangle),
        "--cut",
        "1e300",
    ]
    proc = subprocess.run(
        cmd,
        input=(sequence + "\n").encode("utf-8"),
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=True,
    )
    trajectories = parse_cli_trajectories(proc.stdout.decode("utf-8"))
    if len(trajectories) != args.trajectories_per_sequence:
        raise RuntimeError(
            f"{binary} returned {len(trajectories)} completed trajectories, "
            f"expected {args.trajectories_per_sequence}"
        )
    return trajectories


def run_kinfold_metrics_python(RNA, rows, args):
    rng = random.Random(args.seed)
    seeds = [make_seed_triplet(rng) for _ in rows]
    sequences = [row["sequence"] for row in rows]

    static_mfes = {row["id"]: compute_static_mfe(RNA, row["sequence"]) for row in rows}

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
    for row, seed_triplet, group in zip(rows, seeds, groups):
        case = {
            "case_id": row["id"],
            "case_type": "measured",
            "stem_length": None,
            "sequence": row["sequence"],
            "description": "Measured intrinsic terminator input sequence",
        }
        summary = summarize_case(case, group["trajectories"], static_mfes[row["id"]])
        summary["seed_triplet"] = seed_triplet
        summary["measured_efficiency"] = row["measured_efficiency"]
        case_summaries.append(summary)

    return {"case_summaries": case_summaries, "wall_ns": wall_ns}


def run_kinfold_metrics_cli(RNA, rows, args, binary):
    rng = random.Random(args.seed)
    seeds = [make_seed_triplet(rng) for _ in rows]
    static_mfes = {row["id"]: compute_static_mfe(RNA, row["sequence"]) for row in rows}

    start_ns = time.perf_counter_ns()
    case_summaries = []
    for row, seed_triplet in zip(rows, seeds):
        trajectories = run_kinfold_cli(binary, row["sequence"], seed_triplet, args)
        case = {
            "case_id": row["id"],
            "case_type": "measured",
            "stem_length": None,
            "sequence": row["sequence"],
            "description": "Measured intrinsic terminator input sequence",
        }
        summary = summarize_case(case, trajectories, static_mfes[row["id"]])
        summary["seed_triplet"] = seed_triplet
        summary["measured_efficiency"] = row["measured_efficiency"]
        case_summaries.append(summary)
    wall_ns = time.perf_counter_ns() - start_ns
    return {"case_summaries": case_summaries, "wall_ns": wall_ns}


def feature_value(case, name, time_cap):
    value = case.get(name)
    if value is None:
        return float(time_cap)
    return float(value)


def build_feature_matrix(case_summaries, feature_names, time_cap):
    matrix = []
    target = []
    ids = []
    for case in case_summaries:
        matrix.append([feature_value(case, name, time_cap) for name in feature_names])
        target.append(float(case["measured_efficiency"]))
        ids.append(case["case_id"])
    return np.asarray(matrix, dtype=float), np.asarray(target, dtype=float), ids


def make_cv(n_samples, requested_folds, seed):
    folds = min(requested_folds, n_samples)
    if folds < 2:
        raise ValueError("Need at least 2 samples for cross-validation")
    return KFold(n_splits=folds, shuffle=True, random_state=seed)


def fit_lasso(X, y, feature_names, seed, cv_folds, max_iter):
    outer_cv = make_cv(len(y), cv_folds, seed)
    inner_cv = make_cv(len(y), cv_folds, seed)

    model = Pipeline(
        steps=[
            ("scaler", StandardScaler()),
            ("lasso", LassoCV(cv=inner_cv, random_state=seed, max_iter=max_iter)),
        ]
    )
    model.fit(X, y)

    fitted = model.named_steps["lasso"]
    predictions = model.predict(X)

    cv_model = Pipeline(
        steps=[
            ("scaler", StandardScaler()),
            ("lasso", Lasso(alpha=float(fitted.alpha_), max_iter=max_iter)),
        ]
    )
    cv_predictions = cross_val_predict(cv_model, X, y, cv=outer_cv)

    coefs = fitted.coef_
    selected = [
        {"feature": name, "coefficient": float(coef), "abs_coefficient": float(abs(coef))}
        for name, coef in zip(feature_names, coefs)
        if abs(coef) > 1e-9
    ]
    selected.sort(key=lambda item: item["abs_coefficient"], reverse=True)

    pearson = pearsonr(y, predictions) if len(y) >= 2 else (None, None)
    spearman = spearmanr(y, predictions) if len(y) >= 2 else (None, None)
    cv_pearson = pearsonr(y, cv_predictions) if len(y) >= 2 else (None, None)
    cv_spearman = spearmanr(y, cv_predictions) if len(y) >= 2 else (None, None)

    return {
        "pipeline": model,
        "predictions": predictions,
        "cv_predictions": cv_predictions,
        "selected_features": selected,
        "all_coefficients": [
            {"feature": name, "coefficient": float(coef)}
            for name, coef in zip(feature_names, coefs)
        ],
        "alpha": float(fitted.alpha_),
        "intercept": float(fitted.intercept_),
        "train_r2": float(r2_score(y, predictions)),
        "train_rmse": float(math.sqrt(mean_squared_error(y, predictions))),
        "cv_r2": float(r2_score(y, cv_predictions)),
        "cv_rmse": float(math.sqrt(mean_squared_error(y, cv_predictions))),
        "train_pearson_r": float(pearson.statistic),
        "train_pearson_p": float(pearson.pvalue),
        "train_spearman_rho": float(spearman.statistic),
        "train_spearman_p": float(spearman.pvalue),
        "cv_pearson_r": float(cv_pearson.statistic),
        "cv_pearson_p": float(cv_pearson.pvalue),
        "cv_spearman_rho": float(cv_spearman.statistic),
        "cv_spearman_p": float(cv_spearman.pvalue),
    }


def build_output_rows(case_summaries, feature_names, predictions, cv_predictions, time_cap):
    rows = []
    for case, pred, cv_pred in zip(case_summaries, predictions, cv_predictions):
        row = {
            "case_id": case["case_id"],
            "sequence": case["sequence"],
            "measured_efficiency": case["measured_efficiency"],
            "predicted_efficiency": float(pred),
            "cv_predicted_efficiency": float(cv_pred),
        }
        for feature_name in feature_names:
            row[feature_name] = feature_value(case, feature_name, time_cap)
        rows.append(row)
    return rows


def write_output_csv(path, rows, feature_names):
    fieldnames = [
        "case_id",
        "sequence",
        "measured_efficiency",
        "predicted_efficiency",
        "cv_predicted_efficiency",
        *feature_names,
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    args = parse_args()
    rows = load_rows(args.input_csv, args.id_column, args.sequence_column, args.efficiency_column)

    import RNA

    repo_root = pathlib.Path(__file__).resolve().parents[2]
    binary = args.binary or str(repo_root / "src" / "Kinfold" / "Kinfold")

    if args.backend == "python":
        if not hasattr(RNA, "kinfold_batch"):
            raise RuntimeError("Python backend requested but RNA.kinfold_batch is unavailable")
        run = run_kinfold_metrics_python(RNA, rows, args)
    elif args.backend == "cli":
        run = run_kinfold_metrics_cli(RNA, rows, args, binary)
    else:
        if hasattr(RNA, "kinfold_batch"):
            run = run_kinfold_metrics_python(RNA, rows, args)
        else:
            run = run_kinfold_metrics_cli(RNA, rows, args, binary)

    case_summaries = run["case_summaries"]
    X, y, ids = build_feature_matrix(case_summaries, args.feature_columns, args.time)
    fit = fit_lasso(X, y, args.feature_columns, args.seed, args.cv_folds, args.max_iter)
    output_rows = build_output_rows(case_summaries, args.feature_columns, fit["predictions"], fit["cv_predictions"], args.time)

    payload = {
        "input_csv": str(pathlib.Path(args.input_csv).resolve()),
        "module_path": getattr(RNA, "__file__", None),
        "sequence_count": len(rows),
        "feature_count": len(args.feature_columns),
        "feature_names": args.feature_columns,
        "backend": args.backend,
        "binary": binary,
        "trajectories_per_sequence": args.trajectories_per_sequence,
        "time": args.time,
        "time_step": args.time_step,
        "transcription_elongation_rate": args.transcription_elongation_rate,
        "kinfold_seconds_per_time_unit": args.kinfold_seconds_per_time_unit,
        "max_bubble_width": args.max_bubble_width,
        "temperature": args.temperature,
        "dangle": args.dangle,
        "seed": args.seed,
        "cv_folds": min(args.cv_folds, len(rows)),
        "kinfold_wall_ns": run["wall_ns"],
        "lasso": {
            key: value
            for key, value in fit.items()
            if key not in {"pipeline", "predictions", "cv_predictions"}
        },
        "cases": output_rows,
    }

    pathlib.Path(args.output_json).write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    if args.output_csv:
        write_output_csv(args.output_csv, output_rows, args.feature_columns)


if __name__ == "__main__":
    main()

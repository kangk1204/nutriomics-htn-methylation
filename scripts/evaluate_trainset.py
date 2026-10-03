"""Fresh full-development refits scored on the same participants (resubstitution)."""
from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import warnings
from datetime import datetime, timezone
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.exceptions import ConvergenceWarning
from threadpoolctl import threadpool_limits

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from nutriomics_methylation.evaluation import point_metrics
from nutriomics_methylation.io import sha256, versions, write_json


def record(path):
    return {"path": str(path.resolve()), "sha256": sha256(path), "bytes": path.stat().st_size}


def run(root: Path, output: Path):
    if output.exists():
        raise FileExistsError("Use a new output directory; historical results must be preserved")
    prepared = root / "data/prepared/matrices.npz"
    inputs = [prepared, root / "results/training_config.json", root / "docs/validation/metrics_all.json"]
    inputs += sorted((root / "results").glob("*/final_model.joblib"))
    inputs += sorted((root / "results").glob("*/final_training_audit.json"))
    inputs += sorted((root / "src/nutriomics_methylation").glob("*.py"))
    before = [record(p) for p in inputs]
    config = json.loads((root / "results/training_config.json").read_text())
    if sha256(prepared) != config["prepared_matrix_sha256"]:
        raise ValueError("Prepared input differs from the original frozen training input")
    data = np.load(prepared, allow_pickle=False)
    phenotypes = data["primary_phenotypes"]
    include = phenotypes != "preHT"
    y = (phenotypes[include] == "HTN").astype(np.int8)
    samples = data["primary_samples"][include]
    if len(y) != 88 or int(y.sum()) != 44 or len(np.unique(samples)) != 88:
        raise ValueError("Expected exactly 44 HTN and 44 controls, unique public participants")
    if int(np.sum(~include)) != 44:
        raise ValueError("The 44 preHT participants must remain excluded")
    probes = data["probe_ids"]
    primary = data["primary"][include]
    original_metrics = json.loads((root / "docs/validation/metrics_all.json").read_text())
    output.mkdir(parents=True)
    rows = []
    for space in ("common", "full"):
        indices = np.arange(len(probes)) if space == "full" else pd.Index(probes).get_indexer(data["common_probes"])
        if np.any(indices < 0):
            raise ValueError("Missing common probe identifiers")
        X = primary[:, indices]
        for family in ("ridge", "elasticnet"):
            name = f"{space}_{family}"
            frozen_path = root / "results" / name / "final_model.joblib"
            frozen = joblib.load(frozen_path)
            if not np.array_equal(frozen["training_samples"], samples) or not np.array_equal(frozen["probes"], probes[indices]):
                raise ValueError("Frozen training participants or probe order differ")
            fitted = clone(frozen["model"])
            with threadpool_limits(limits=1), warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always", ConvergenceWarning)
                fitted.fit(X, y)
                probability = fitted.predict_proba(X)[:, 1]
                prior_probability = frozen["model"].predict_proba(X)[:, 1]
            metrics = point_metrics(y, probability)
            metrics.update(accuracy=float(np.mean((probability >= 0.5) == y)),
                           evaluation="training_resubstitution", independent_evaluation=False)
            pd.DataFrame({"sample_id": samples, "y": y, "probability": probability,
                          "used_for_training": True, "used_for_evaluation": True}).to_csv(
                output / f"{name}_training_predictions.tsv", sep="\t", index=False)
            tuning_path = root / "results" / name / "final_training_audit.json"
            tuning = json.loads(tuning_path.read_text())
            details = {"model": name, "metrics": metrics, "parameters": tuning["best_params"],
                       "parameters_source": "Original inner-CV selected final refit parameters, fixed before this diagnostic",
                       "selected_probes": len(fitted["probes"].selected_indices_),
                       "available_probes": int(X.shape[1]),
                       "fresh_fit_n": int(fitted["probes"].n_fit_samples_),
                       "refit_convergence_warnings": sum(isinstance(w.message, ConvergenceWarning) for w in caught),
                       "original_refit_max_probability_difference": float(np.max(np.abs(probability-prior_probability))),
                       "nested5_auroc": original_metrics[name + "/nested5"]["auroc"],
                       "loocv_auroc": original_metrics[name + "/loocv"]["auroc"],
                       "external_auroc": original_metrics.get(name + "/external", {}).get("auroc")}
            write_json(output / f"{name}_metrics.json", details)
            rows.append(details)
            print(json.dumps({"model": name, "train_auroc": metrics["auroc"],
                              "nested5_auroc": details["nested5_auroc"], "accuracy": metrics["accuracy"]}), flush=True)
    after = [record(p) for p in inputs]
    if before != after:
        raise ValueError("Historical inputs, code or primary results changed during the experiment")
    table = [{"model": r["model"], "training_auroc": r["metrics"]["auroc"],
              "nested5_auroc": r["nested5_auroc"], "loocv_auroc": r["loocv_auroc"],
              "external_auroc": r["external_auroc"], "training_accuracy": r["metrics"]["accuracy"],
              "training_brier": r["metrics"]["brier"], "selected_probes": r["selected_probes"]} for r in rows]
    pd.DataFrame(table).to_csv(output / "comparison.tsv", sep="\t", index=False)
    try:
        revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True, stderr=subprocess.DEVNULL).strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        revision = None
    write_json(output / "receipt.json", {
        "generated_utc": datetime.now(timezone.utc).isoformat(), "host": platform.node(),
        "versions": versions(), "code_revision_at_run": revision,
        "script": record(Path(__file__)), "source_files": before,
        "historical_sources_unchanged": True, "participants": 88, "htn": 44, "controls": 44,
        "preHT_excluded": 44, "random_seed": config["seed"],
        "fresh_full_development_refits": 4, "same_training_and_evaluation_participants": True,
        "purpose": "Requested diagnostic of optimism from fitting and scoring the same cohort",
        "clinical_or_independent_performance_claim": False,
        "historical_auc_088_protocol_confirmed": False,
        "limits": "Does not reconstruct unavailable historical 131/678 participants or prove the historical evaluation protocol. No generalization or clinical accuracy claim.",
        "results": rows,
        "output_files": [record(p) for p in sorted(output.iterdir()) if p.is_file()],
    })
    lines = ["# 고혈압 메틸화 학습자료 재평가", "",
             "GSE193795 고혈압44명·대조44명을 모두 학습한 뒤 같은88명을 평가했다. 전고혈압44명은 제외했다.",
             "기존 최종 모델의 내부 교차검증 선택 파라미터를 고정해 새로4회 적합했다. 학습 AUROC는 resubstitution 적합도이며 독립 검증 성능이 아니다.", "",
             "|모델|학습자료 AUROC|nested5 AUROC|nested LOO AUROC|외부 AUROC|", "|---|---:|---:|---:|---:|"]
    for r in rows:
        ext = f"{r['external_auroc']:.3f}" if r["external_auroc"] is not None else "해당 없음"
        lines.append(f"|{r['model']}|{r['metrics']['auroc']:.3f}|{r['nested5_auroc']:.3f}|{r['loocv_auroc']:.3f}|{ext}|")
    lines += ["", "학습자료에서 높은 값이 나와도 기존AUC0.88의 평가 방식은 원본 모델·분할·예측값 없이는 확정할 수 없다. 새 학습 적합도로 과거 실적이나 독립 평가 결과를 대체하지 않는다.", ""]
    (output / "SUMMARY_KO.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.root.resolve(), args.output.resolve())

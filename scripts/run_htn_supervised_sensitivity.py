"""Exploratory F-statistic HTN selector, compared with preserved variance models.

The small external cohort has already been evaluated by the original models;
this additional analysis is explicitly post-hoc, not an untouched new test.
"""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import os
from pathlib import Path
import sys
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(variable, "1")
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from nutriomics_methylation.chronic import nested, tune
from nutriomics_methylation.evaluation import performance, point_metrics, match_groups
from nutriomics_methylation.io import sha256, versions, write_json


def run(root, output):
    if output.exists():
        raise FileExistsError("Use new output paths")
    prepared = root / "data/prepared/matrices.npz"
    digest = sha256(prepared)
    output.mkdir(parents=True)
    seeds = [20261001, 20261018, 20261101]
    write_json(output / "protocol.json", {
        "posthoc": True, "seeds": seeds, "primary_n": 88, "htn_n": 44, "control_n": 44,
        "preHT_excluded_from_learning": 44, "selector": "source-fold-local F statistic",
        "nested": "5 outer / 3 inner", "k": [10, 50, 200], "C": [.1, 1.],
        "spaces": ["common", "full"], "threshold": .5, "prepared_sha256": digest,
        "external_already_seen_in_previous_analyses": True, "external_tuning": False,
        "age_sex_adjustment": False, "clinical_validation": False})
    data = np.load(prepared, allow_pickle=False)
    include = data["primary_phenotypes"] != "preHT"
    y = (data["primary_phenotypes"][include] == "HTN").astype(int)
    samples = data["primary_samples"][include]
    if np.bincount(y).tolist() != [44, 44] or (~include).sum() != 44:
        raise ValueError("Frozen HTN phenotype denominators changed")
    external_meta = pd.read_csv(root / "data/prepared/external_metadata.tsv", sep="\t", index_col=0)
    if external_meta.index.tolist() != data["external_samples"].tolist():
        raise ValueError("External participant metadata differs")
    groups = match_groups(external_meta)
    rows = []
    with threadpool_limits(limits=1):
        for space in ("common", "full"):
            dest = output / space
            dest.mkdir()
            probes = data["common_probes"] if space == "common" else data["probe_ids"]
            indices = pd.Index(data["probe_ids"]).get_indexer(probes)
            if (indices < 0).any():
                raise ValueError("Probe alignment failed")
            X = data["primary"][include][:, indices].astype(np.float32)
            probabilities = []
            for seed in seeds:
                pred, fold, audits = nested(X, y, samples, seed, "supervised")
                pd.DataFrame({"sample_id": samples, "y": y, "probability": pred, "fold": fold}).to_csv(dest / f"{seed}_oof.tsv", sep="\t", index=False)
                write_json(dest / f"{seed}_fold_audit.json", audits)
                probabilities.append(pred)
                print(f"space={space} seed={seed} AUROC={point_metrics(y, pred)['auroc']:.6f}", flush=True)
            pred = np.mean(probabilities, axis=0)
            pd.DataFrame({"sample_id": samples, "y": y, "probability": pred}).to_csv(dest / "repeated_oof.tsv", sep="\t", index=False)
            fitted, audit = tune(X, y, seeds[0], "supervised")
            write_json(dest / "source_only_tuning.json", audit)
            pd.DataFrame({"feature_id": probes[fitted["probes"].selected_indices_], "standardized_coefficient": fitted["model"].coef_[0]}).to_csv(dest / "refit_features.tsv", sep="\t", index=False)
            result = {"space": space, "source": "GSE193795", "target": "GSE42774" if space == "common" else None,
                      "source_metrics": performance(y, pred, seed=seeds[0]), "replication_metrics": None,
                      "training_resubstitution": point_metrics(y, fitted.predict_proba(X)[:, 1]),
                      "selected_k": audit["k"], "selected_C": audit["C"], "posthoc": True,
                      "target_fitting_tuning_calibration_or_threshold_selection": False}
            if space == "common":
                ext = fitted.predict_proba(data["external"])[:, 1]
                pd.DataFrame({"sample_id": data["external_samples"], "y": data["external_y"], "probability": ext}).to_csv(dest / "external.tsv", sep="\t", index=False)
                result["replication_metrics"] = performance(data["external_y"], ext, seed=seeds[0], groups=groups)
            write_json(dest / "summary.json", result)
            rows.append(result)
    if sha256(prepared) != digest:
        raise ValueError("Original prepared data changed")
    write_json(output / "summary.json", {"analyses": rows, "original_matrix_unchanged": True,
               "versions": versions(), "script_sha256": sha256(Path(__file__)), "generated_utc": datetime.now(timezone.utc).isoformat()})
    print("HTN supervised diagnostic completed", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    run(args.root, args.output)

"""Post-hoc SNP/probe-annotation robustness; original frozen result is preserved.

The entire native 450K platform is screened, rather than deleting only the
eight flagged features in the original final fit. Selection, imputation,
scaling and tuning remain within source folds. Replication is never tuned.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
from pathlib import Path
import sys
from datetime import datetime, timezone

for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(variable, "1")
import joblib
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from nutriomics_methylation.chronic import nested, tune
from nutriomics_methylation.evaluation import performance, point_metrics
from nutriomics_methylation.io import sha256, versions, write_json

SEEDS = (20261003, 20261020, 20261103)


def flagged_platform(path):
    flagged, seen, active, header = set(), set(), False, None
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        for line in stream:
            if line.startswith("!platform_table_begin"):
                active = True
                continue
            if line.startswith("!platform_table_end"):
                active = False
                continue
            if not active:
                continue
            parts = next(csv.reader([line], delimiter="\t", quoting=csv.QUOTE_NONE))
            if header is None:
                header = parts
                for column in ("ID", "Probe_SNPs", "Probe_SNPs_10"):
                    if column not in header:
                        raise ValueError("Native probe/SNP columns are missing")
                continue
            if len(parts) != len(header):
                raise ValueError("Native annotation width mismatch")
            row = dict(zip(header, parts))
            probe = row["ID"]
            if probe in seen:
                raise ValueError("Duplicate platform probe ID")
            seen.add(probe)
            if row["Probe_SNPs"].strip() or row["Probe_SNPs_10"].strip():
                flagged.add(probe)
    if not seen or active:
        raise ValueError("Incomplete native platform table")
    return flagged, seen


def run(root, output):
    if output.exists():
        raise FileExistsError("Preserve earlier results; use a new output directory")
    prepared = root / "data/chronic_20261003/prepared_v3"
    annotation = root / "data/chronic_20261003/obesity_annotation_v1/GPL13534_self.soft.gz"
    original_model = root / "models/chronic_20261003_obesity_discovery_supervised.joblib"
    originals = [prepared / (name + ".npz") for name in ("obesity_discovery", "obesity_replication")]
    originals += [annotation, original_model]
    before = {str(p.relative_to(root)): sha256(p) for p in originals}
    output.mkdir(parents=True)
    write_json(output / "protocol.json", {
        "analysis": "post-hoc annotation-based robustness", "locked_before_this_run_results": True,
        "original_replication_results_already_known": True, "seeds": list(SEEDS),
        "exclusion_rule": "Any nonempty native Probe_SNPs or Probe_SNPs_10 on the entire GPL13534 platform; no MAF threshold claimed",
        "nested": "5 outer / 3 inner folds, source-only learned filtering/imputation/F-statistic ranking/scaling/tuning",
        "grid": {"k": [10, 50, 200], "C": [0.1, 1.0]}, "threshold": 0.5,
        "replication_training_or_selection": False, "clinical_validation": False,
        "raw_IDAT_reprocessing_or_covariate_adjustment": False,
        "all_inputs_sha256": before})
    flagged, annotated = flagged_platform(annotation)
    source = dict(np.load(originals[0], allow_pickle=False))
    target = dict(np.load(originals[1], allow_pickle=False))
    probes = np.array(sorted(set(source["probes"]) & set(target["probes"])))
    lookup_a = {p: i for i, p in enumerate(source["probes"])}
    lookup_b = {p: i for i, p in enumerate(target["probes"])}
    X = source["X"][:, [lookup_a[p] for p in probes]]
    Z = target["X"][:, [lookup_b[p] for p in probes]]
    y, z, samples = source["y"], target["y"], source["samples"]
    if np.bincount(y).tolist() != [24, 24] or np.bincount(z).tolist() != [24, 24]:
        raise ValueError("Fixed cohort denominators changed")
    if set(samples) & set(target["samples"]):
        raise ValueError("Source-target accession overlap")
    keep = np.array([p not in flagged for p in probes])
    removed = probes[~keep]
    pd.DataFrame({"feature_id": removed}).to_csv(output / "excluded_native_flagged_probes.tsv", sep="\t", index=False)
    summary = {"platform_rows": len(annotated), "platform_flagged": len(flagged),
               "aligned_union_probes": len(probes), "removed_aligned_probes": len(removed),
               "remaining_aligned_probes": int(keep.sum()),
               "source_usable_90pct_before": int((np.mean(np.isfinite(X), axis=0) >= .9).sum()),
               "source_usable_90pct_after": int((np.mean(np.isfinite(X[:, keep]), axis=0) >= .9).sum())}
    Xa, Za, retained = X[:, keep], Z[:, keep], probes[keep]
    with threadpool_limits(limits=1):
        probabilities = []
        for seed in SEEDS:
            pred, fold, audit = nested(Xa, y, samples, seed, "supervised")
            probabilities.append(pred)
            pd.DataFrame({"sample_id": samples, "y": y, "probability": pred, "fold": fold}).to_csv(output / f"{seed}_oof.tsv", sep="\t", index=False)
            write_json(output / f"{seed}_fold_audit.json", audit)
            print(json.dumps({"seed": seed, "source_auroc": point_metrics(y, pred)["auroc"]}), flush=True)
        pred = np.mean(probabilities, axis=0)
        fitted, audit = tune(Xa, y, SEEDS[0], "supervised")
        ext = fitted.predict_proba(Za)[:, 1]
        pd.DataFrame({"sample_id": samples, "y": y, "probability": pred}).to_csv(output / "repeated_oof.tsv", sep="\t", index=False)
        pd.DataFrame({"sample_id": target["samples"], "y": z, "probability": ext}).to_csv(output / "external.tsv", sep="\t", index=False)
        selected = retained[fitted["probes"].selected_indices_]
        if set(selected) & flagged:
            raise ValueError("Flagged probe survived exclusion")
        pd.DataFrame({"feature_id": selected, "standardized_coefficient": fitted["model"].coef_[0]}).to_csv(output / "refit_features.tsv", sep="\t", index=False)
        write_json(output / "source_only_tuning.json", audit)
        # Models remain in ignored models/; derived metrics are exportable.
        model_file = root / "models" / (output.name + "_source_only.joblib")
        if model_file.exists():
            raise FileExistsError("Sensitivity model artifact already exists")
        joblib.dump({"model": fitted, "probe_ids": retained, "training_samples": samples}, model_file)
        summary.update(source_metrics=performance(y, pred, seed=SEEDS[0]),
                       replication_metrics=performance(z, ext, seed=SEEDS[0]),
                       training_resubstitution=point_metrics(y, fitted.predict_proba(Xa)[:, 1]),
                       selected_probes=len(selected), selected_k=audit["k"], selected_C=audit["C"],
                       model_sha256=sha256(model_file), model_path=str(model_file))
    if before != {str(p.relative_to(root)): sha256(p) for p in originals}:
        raise ValueError("Original input or frozen model changed")
    summary.update(original_inputs_and_model_unchanged=True, posthoc=True,
                   replication_used_for_tuning=False, versions=versions(),
                   script_sha256=sha256(Path(__file__)), generated_utc=datetime.now(timezone.utc).isoformat())
    write_json(output / "summary.json", summary)
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    run(args.root, args.output)

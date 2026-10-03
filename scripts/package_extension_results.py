"""Package executed extension results without promoting prior results to new runs."""
import argparse
import hashlib
import json
from pathlib import Path
from datetime import datetime, timezone
import shutil
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def record(path):
    return {"path": str(path.relative_to(ROOT)), "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def metric_row(task, model, space, cohort, role, metrics, receipt):
    result = {"task": task, "model": model, "feature_space": space, "cohort": cohort, "evaluation_role": role,
              "n": metrics["n"], "n_positive": metrics["n_positive"],
              "receipt_path": str(receipt.relative_to(ROOT)), "receipt_sha256": record(receipt)["sha256"]}
    for key in ["auroc", "average_precision", "auprc_trapezoid", "brier", "sensitivity", "specificity",
                "threshold", "calibration_slope_descriptive", "calibration_intercept_descriptive"]:
        result[key] = metrics.get(key)
    result["auroc_ci_low"], result["auroc_ci_high"] = metrics["ci95"]["auroc"]
    result["interval_method"] = metrics["interval_method"]["discrimination_brier"]
    result["interval_scope"] = metrics["interval_method"]["scope"]
    return result


def main():
    parser = argparse.ArgumentParser(); parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--include-endocrine", action="store_true", help="Add all16 audited rows in a distinct combined package")
    args = parser.parse_args()
    out = args.output.resolve(); out.mkdir(parents=True, exist_ok=True)
    if args.include_endocrine and out == (ROOT / "docs/phase1_extension_receipts_20261003").resolve():
        raise ValueError("The original14 output path is reserved; use a separate directory for the combined30 package")
    rows = []; inputs = []
    obesity_path = ROOT / "docs/obesity_replay_20261003_v1/summary.json"
    obesity = json.loads(obesity_path.read_text()); inputs.append(obesity_path)
    roles = {"source_nested": "independent source nested-CV replay, not a fixed 50-probe CV model",
             "source_resubstitution": "training-set diagnostic; no generalization claim",
             "old_epic_replication": "already-seen independent computational replay; upstream jointly processed matrix",
             "GSE67024_seen_replay": "already-seen separate-cohort computational replay; no target fit"}
    for value in obesity["results"]:
        rows.append(metric_row("Obese versus lean adipocytes", "source logistic grid" if value["cohort"] == "source_nested" else "frozen obesity50 logistic",
                               "source-selected", value["cohort"], roles[value["cohort"]], value["metrics"], obesity_path))
    twins_path = ROOT / "docs/twins_transfer_20261003_v2/summary.json"
    twins = json.loads(twins_path.read_text()); inputs.append(twins_path)
    rows.append(metric_row("Within-family higher versus lower BMI", "frozen obesity50 logistic", "source-selected",
                           "GSE68336 bulk SAT", "cross-tissue endpoint transfer; 24 discordant MZ families; not clinical obesity",
                           twins["primary_relative_BMI_transfer_metrics"], twins_path))
    for space in ["common", "full"]:
        for model in ["ridge", "elasticnet", "linear_svm"]:
            path = ROOT / f"docs/htn_grid_20261003_v1/{space}/{model}_summary.json"
            item = json.loads(path.read_text()); inputs.append(path)
            rows.append(metric_row("HTN versus control", model, space, "GSE193795", "source-only nested CV; source-OOF sigmoid calibration",
                                   item["source_nested"], path))
            if space == "common":
                rows.append(metric_row("HTN versus control", model, space, "GSE42774", "already-seen matched-pair target diagnostic; no target tuning",
                                       item["seen_external_diagnostic"], path))
    endocrine_dir = ROOT / "docs/endocrine_extension_20261003_v1"
    endocrine_summary = endocrine_dir / "summary.json"
    endocrine_audit = endocrine_dir / "audit.json"
    endocrine_evaluated = False
    endocrine_raw_receipt = ROOT / "docs/endocrine_acquisition_complete_20261003_v1/collection_receipt.json"
    endocrine_raw_complete = endocrine_raw_receipt.exists() and json.loads(endocrine_raw_receipt.read_text())["whole_raw_verified_complete"]
    if endocrine_raw_complete:
        inputs.append(endocrine_raw_receipt)
    if endocrine_summary.exists() and endocrine_audit.exists():
        summary = json.loads(endocrine_summary.read_text())
        audit = json.loads(endocrine_audit.read_text())
        expected = {(endpoint, family) for endpoint in ["CS_vs_PHT", "PA_vs_PHT", "PPGL_vs_PHT", "EHT_vs_PHT"]
                    for family in ["ridge", "elasticnet", "linear_svm", "age_sex_only"]}
        actual = {(r["endpoint"], r["family"]) for r in summary["analyses"]}
        if not audit["all_16_results_passed"] or audit["summary_sha256"] != record(endocrine_summary)["sha256"]:
            raise ValueError("Endocrine summary has not passed the current independent audit")
        if actual != expected or len(summary["analyses"]) != 16 or not summary["all_four_endpoints_reported"] or summary["healthy_control_count"] != 0:
            raise ValueError("All declared endocrine endpoints/families and control definition required")
        labels = {"CS_vs_PHT": "CS versus primary HTN", "PA_vs_PHT": "PA versus primary HTN",
                  "PPGL_vs_PHT": "PPGL versus primary HTN", "EHT_vs_PHT": "Combined endocrine CS+PA+PPGL versus primary HTN"}
        for value in summary["analyses"]:
            family = value["family"]
            if args.include_endocrine:
                rows.append(metric_row(labels[value["endpoint"]], family,
                                       "age/sex comparator" if family == "age_sex_only" else "repository raw-beta CpGs",
                                       "S-BSST831", "same-study nested center-held-out evaluation; no healthy controls; center/subtype confounding",
                                       value["nested_center_metrics"], endocrine_summary))
        inputs.extend([endocrine_summary, endocrine_audit])
        endocrine_evaluated = True
    pd.DataFrame(rows).to_csv(out / "all_executed_extension_evaluations.tsv", sep="\t", index=False)
    prior = out / "prior_comparisons_not_newly_reexecuted"; prior.mkdir(exist_ok=True)
    for source in [ROOT / "docs/chronic_20261003/screen_v1/comparison.tsv",
                   ROOT / "docs/chronic_20261003/methylation_transfer_v3/comparison.tsv",
                   ROOT / "docs/trainset_20261003/comparison.tsv"]:
        inputs.append(source); shutil.copy2(source, prior / f"{source.parent.name}_comparison.tsv")
    coverage = [{"protocol": "obesity immutable source model replay", "state": "completed", "evidence": "docs/obesity_replay_20261003_v1/summary.json"},
                {"protocol": "HTN nested ridge/elastic-net/linear-SVM", "state": "completed; AUC >0.8 not met", "evidence": "docs/htn_grid_20261003_v1/summary.json"},
                {"protocol": "GSE68336 family-preserved frozen-model transfer", "state": "completed; cross-tissue diagnostic", "evidence": "docs/twins_transfer_20261003_v2/summary.json"},
                {"protocol": "S-BSST831 raw-beta endocrine contrasts", "state": "completed; all16 independently audited; no healthy controls" if endocrine_evaluated else
                 ("raw complete; evaluation or independent audit in progress" if endocrine_raw_complete else "acquisition pending; no model results"),
                 "evidence": "docs/endocrine_extension_20261003_v1/audit.json" if endocrine_evaluated else
                 ("docs/endocrine_acquisition_complete_20261003_v1/collection_receipt.json" if endocrine_raw_complete else "scripts/run_endocrine_extension.py")},
                {"protocol": "GSE314180", "state": "mechanistic only; no verified healthy controls", "evidence": "docs/extension_eligibility_20261003_v1/eligibility.json"},
                {"protocol": "GSE262817", "state": "no individual numeric matrix acquired; no DMR classifier", "evidence": "docs/extension_eligibility_20261003_v1/eligibility.json"}]
    pd.DataFrame(coverage).to_csv(out / "protocol_status.tsv", sep="\t", index=False)
    receipt = {"generated_at": datetime.now(timezone.utc).isoformat(), "evaluations": len(rows), "input_records": [record(p) for p in inputs],
               "obesity_original_predictions_reproduced": obesity["all_existing_predictions_reproduced"],
               "obesity_model_unchanged": obesity["frozen_model_unchanged"], "obesity_model_sha256": twins["model_sha256"],
               "HTN_target_0_8_met": False, "endocrine_model_evaluated": endocrine_evaluated,
               "endocrine_raw_complete": bool(endocrine_raw_complete),
               "endocrine_rows_in_this_table": 16 if args.include_endocrine and endocrine_evaluated else 0,
               "clinical_GSE68336_obesity_claim": False, "prior_comparison_tables_are_new_executions": False,
               "regenerate": "python scripts/package_extension_results.py --output " + str(out.relative_to(ROOT)) +
               (" --include-endocrine" if args.include_endocrine else "")}
    (out / "package_receipt.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    print(json.dumps({"evaluations": len(rows), "output": str(out)}))


if __name__ == "__main__": main()

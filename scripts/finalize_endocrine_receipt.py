"""Produce a small report receipt from completed, independently audited outputs."""
from __future__ import annotations

import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def record(path: Path):
    return {"path": str(path.relative_to(ROOT)).replace("\\", "/"),
            "bytes": path.stat().st_size,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def write(path: Path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def main():
    final = ROOT / "docs/endocrine_extension_20261003_v1"
    original = ROOT / "docs/phase1_extension_receipts_20261003"
    raw = ROOT / "docs/endocrine_acquisition_complete_20261003_v1"
    summary = load(final / "summary.json")
    audit = load(final / "audit.json")
    acquisition = load(raw / "range_acquisition_status.json")
    collected = load(raw / "collection_receipt.json")
    with (original / "all_executed_extension_evaluations.tsv").open(encoding="utf-8", newline="") as handle:
        original_rows = list(csv.DictReader(handle, delimiter="\t"))
    with (final / "comparison.tsv").open(encoding="utf-8", newline="") as handle:
        final_rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(original_rows) == 14, "Original report table must remain14 rows"
    assert len(final_rows) == len(summary["analyses"]) == len(audit["analyses"]) == 16
    assert audit["all_16_results_passed"] and audit["no_model_fit_or_selection"]
    assert audit["summary_sha256"] == record(final / "summary.json")["sha256"]
    assert summary["all_four_endpoints_reported"] and audit["healthy_control_count"] == 0
    assert acquisition["state"] == "complete" and collected["whole_raw_verified_complete"]
    assert acquisition["sha256"] == audit["whole_raw_sha256"] == collected["whole_raw_sha256"]
    assert acquisition["verified_bytes"] == acquisition["official_bytes"] == 3381823154
    methylation = [v for v in summary["analyses"] if v["family"] != "age_sex_only"]
    best = max(methylation, key=lambda v: v["nested_center_metrics"]["auroc"])
    now = datetime.now(timezone.utc).isoformat()
    receipt = {
        "observed_utc": now,
        "state": "complete_and_independently_audited",
        "whole_raw_verified_complete": True,
        "partial_only": False,
        "not_analyzed": False,
        "analysis_results_verified": True,
        "verified_bytes": acquisition["verified_bytes"],
        "official_bytes": acquisition["official_bytes"],
        "verified_percent": 100.0,
        "raw_complete_utc": acquisition["completed_utc"],
        "independent_audit_utc": audit["audited_utc"],
        "whole_raw_sha256": acquisition["sha256"],
        "whole_raw_crc32": acquisition["crc32"],
        "locally_generated_checksums_not_official_published_hashes": True,
        "verified_ranges": 404,
        "unique_source_participants": 255,
        "source_cpgs": 731635,
        "healthy_control_count": 0,
        "shared_primary_htn_controls": 42,
        "completed_methylation_evaluations": 12,
        "completed_age_sex_comparators": 4,
        "all_16_results_passed": True,
        "original_report_table_rows": 14,
        "original_report_table_preserved": True,
        "endocrine_table_rows": 16,
        "best_methylation_result": {
            "endpoint": best["endpoint"], "family": best["family"],
            "n": best["nested_center_metrics"]["n"],
            "auroc": best["nested_center_metrics"]["auroc"],
            "auroc_ci95": best["nested_center_metrics"]["ci95"]["auroc"],
            "brier": best["nested_center_metrics"]["brier"],
        },
        "methylation_auroc_above_0_8_met": any(v["nested_center_metrics"]["auroc"] > .8 for v in methylation),
        "independent_external_study": False,
        "normal_vs_general_hypertension_claim": False,
        "age_sex_adjusted_methylation_claim": False,
        "model_selection_or_fit_by_auditor": False,
        "limitations": [
            "All255 participants have hypertension; no healthy-control classification.",
            "EHT_vs_PHT means combined endocrine CS+PA+PPGL versus primary hypertension.",
            "Same42 primary-hypertension controls recur across contrasts; do not sum participant counts.",
            "Recruitment center and subtype are strongly confounded.",
            "Center-cluster bootstrap conditions on fixed OOF predictions, without retraining or split uncertainty.",
            "Age/sex-only comparator does not adjust methylation models for age and sex.",
            "Source-provided raw-beta export does not permit independent IDAT/detection-P processing.",
        ],
        "files": [record(p) for p in [final / "comparison.tsv", final / "summary.json", final / "audit.json",
                  final / "protocol.json", final / "preparation.json", raw / "collection_receipt.json",
                  raw / "range_acquisition_status.json", original / "all_executed_extension_evaluations.tsv"]],
        "regenerate": "python scripts/finalize_endocrine_receipt.py",
        "receipt_script": record(Path(__file__).resolve()),
    }
    report = original / "endocrine_report_state_20261003_v3.json"
    write(report, receipt)
    status_path = original / "endocrine_supervision_status.json"
    status = load(status_path)
    status.update({"observed_utc": now, "state": receipt["state"], "supervisor_pid": 1350045,
                   "download_pid": 1350046, "original_processes_finished": True,
                   "latest_source_snapshot": "../endocrine_acquisition_complete_20261003_v1/range_acquisition_status.json",
                   "latest_bounded_report": report.name, "analysis_results_verified": True,
                   "completed_evaluations": 16, "next_action": "Integrate separate16 comparison and audited negative results in report"})
    write(status_path, status)
    print(json.dumps({"state": receipt["state"], "receipt": record(report), "original_rows": 14,
                      "endocrine_rows": 16, "best_methylation_auroc": receipt["best_methylation_result"]["auroc"]}))


if __name__ == "__main__":
    main()

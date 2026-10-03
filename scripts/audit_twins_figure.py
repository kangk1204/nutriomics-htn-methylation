"""Independently verify central twin figure quantities; this is not a causal audit."""
import argparse
import json
from pathlib import Path
import shutil
import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control, t

ROOT = Path(__file__).resolve().parents[1]


def main():
    parser = argparse.ArgumentParser(); parser.add_argument("--figure", type=Path, required=True); args = parser.parse_args()
    output = args.figure.resolve()
    table = pd.read_csv(output / "annotated_all50_cpg.tsv", sep="\t")
    summary = json.loads((output / "inputs/twin_summary.json").read_text())
    predictions = pd.read_csv(output / "inputs/twin_predictions.tsv", sep="\t")
    selected = predictions.loc[predictions.discordant]
    assert len(predictions) == 70 and predictions.family_id.nunique() == 35
    assert len(selected) == 48 and selected.family_id.nunique() == 24
    assert selected.groupby("family_id").size().eq(2).all()
    assert selected.groupby("family_id").bmi_group.nunique().eq(2).all()
    truth = selected.bmi_group.eq("high").to_numpy()
    probabilities = selected.probability.to_numpy()
    ranks = pd.Series(probabilities).rank(method="average").to_numpy()
    positives = truth.sum(); negatives = (~truth).sum()
    auc = (ranks[truth].sum() - positives * (positives + 1) / 2) / (positives * negatives)
    brier = np.mean((truth.astype(float) - probabilities) ** 2)
    expected = summary["primary_relative_BMI_transfer_metrics"]
    assert np.isclose(auc, expected["auroc"], atol=1e-14)
    assert np.isclose(brier, expected["brier"], atol=1e-14)
    recalculated_q = false_discovery_control(table.p_value.to_numpy(), method="bh")
    np.testing.assert_allclose(recalculated_q, table.q_value_selected_50_family, atol=1e-14, rtol=1e-12)
    halfwidth = t.ppf(.975, 23) * table.twin_standard_error_pp
    np.testing.assert_allclose(table.twin_delta_beta_pp-halfwidth, table.twin_ci_low_pp, atol=1e-12)
    np.testing.assert_allclose(table.twin_delta_beta_pp+halfwidth, table.twin_ci_high_pp, atol=1e-12)
    directions = int((np.sign(table.twin_delta_beta_pp) == np.sign(table.discovery_delta_beta_pp)).sum())
    hits = int((recalculated_q < .05).sum())
    assert len(table) == 50 and table.feature_id.nunique() == 50 and directions == 48 and hits == 4
    result = {"scope": "Source-table consistency, independent metric/BH/paired-CI recalculation; no causal expression claim",
              "total_people": 70, "total_families": 35, "discordant_people": 48, "discordant_families": 24,
              "auroc_rank_sum_recalculated": float(auc), "brier_recalculated": float(brier),
              "directions_agree": directions, "selected50_BH_hits": hits,
              "paired_t_CI_df": 23, "all_50_loci_retained": True,
              "native_gene_unassigned": int(table.gene_annotation_display.eq("Unassigned").sum()),
              "EPIC_effect_unmeasured_retained": int(table.replication_delta_beta_pp.isna().sum()),
              "checks_passed": True}
    (output / "scientific_source_checks.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    if Path(__file__).resolve() != (output / "audit_twins_figure.py").resolve():
        shutil.copy2(__file__, output / "audit_twins_figure.py")
    print(json.dumps(result))


if __name__ == "__main__": main()

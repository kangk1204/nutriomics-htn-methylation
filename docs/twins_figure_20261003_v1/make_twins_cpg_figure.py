"""Render the predeclared obesity50 twin gene-context figure from frozen receipts.

Run from repository root:
python scripts/make_twins_cpg_figure.py --output docs/twins_figure_20261003_v1
Only supplied derived measurements are plotted. No fitting or inference is done.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import platform
import shutil
from datetime import datetime, timezone

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def assemble(annotation, twin, obesity_summary, twin_summary):
    """One-to-one join and explicit units; reject missing or duplicated loci."""
    if len(annotation) != 50 or len(twin) != 50:
        raise ValueError("Exactly all 50 source-selected CpGs must be supplied")
    if annotation.feature_id.duplicated().any() or twin.feature_id.duplicated().any():
        raise ValueError("Duplicated CpG annotation or twin measurement")
    extra = twin[["feature_id", "paired_families", "delta_beta_high_minus_low", "standard_error",
                  "confidence_low", "confidence_high", "p_value", "q_value_selected_50_family"]]
    table = annotation.merge(extra, on="feature_id", validate="one_to_one", how="outer", indicator=True)
    if not table._merge.eq("both").all():
        raise ValueError("CpG identities do not match")
    table = table.drop(columns="_merge").sort_values("source_coefficient_rank")
    table["gene_annotation_display"] = table.gene_symbols_native.fillna("").replace("", "Unassigned")
    table["gene_relation_scope"] = "Native hg19 annotation association; no expression causality"
    for source, dest in [("delta_beta_high_minus_low", "twin_delta_beta_pp"),
                         ("standard_error", "twin_standard_error_pp"),
                         ("confidence_low", "twin_ci_low_pp"), ("confidence_high", "twin_ci_high_pp")]:
        table[dest] = pd.to_numeric(table[source]) * 100
    table["twin_direction_matches_discovery"] = np.sign(table.twin_delta_beta_pp) == np.sign(table.discovery_delta_beta_pp)
    table["twin_q_below_005_selected50"] = table.q_value_selected_50_family < .05
    if not table.paired_families.eq(24).all():
        raise ValueError("The inferential unit must remain 24 discordant families")
    if not np.isfinite(table[["discovery_delta_beta_pp", "twin_delta_beta_pp"]]).all().all():
        raise ValueError("Cannot silently omit nonfinite CpG effects")
    table["EPIC_effect_measured"] = np.isfinite(table.replication_delta_beta_pp)
    if not twin_summary["no_target_fitting_or_calibration"] or not twin_summary["frozen_model_unchanged"]:
        raise ValueError("Fixed-model receipt does not pass")
    aurocs = []
    labels = {"source_nested": "Nested source: adipocytes", "old_epic_replication": "EPIC replay: adipocytes",
              "GSE67024_seen_replay": "GSE67024 replay: adipocytes"}
    for result in obesity_summary["results"]:
        if result["cohort"] in labels:
            metric = result["metrics"]
            aurocs.append({"cohort": result["cohort"], "label": labels[result["cohort"]],
                           "n": metric["n"], "auroc": metric["auroc"],
                           "ci_low": metric["ci95"]["auroc"][0], "ci_high": metric["ci95"]["auroc"][1],
                           "interval_method": metric["interval_method"]["discrimination_brier"],
                           "clinical_obesity_endpoint": True})
    metric = twin_summary["primary_relative_BMI_transfer_metrics"]
    aurocs.append({"cohort": "GSE68336_bulk_SAT_transfer", "label": "Twin transfer: bulk SAT", "n": metric["n"],
                   "auroc": metric["auroc"], "ci_low": metric["ci95"]["auroc"][0], "ci_high": metric["ci95"]["auroc"][1],
                   "interval_method": metric["interval_method"]["discrimination_brier"], "clinical_obesity_endpoint": False})
    if len(aurocs) != 4:
        raise ValueError("Missing AUROC receipts")
    return table, pd.DataFrame(aurocs)


def make_figure(table, aurocs, output):
    # Nature specifications: editable text, restrained spines, 183-mm width, >=5-pt labels.
    plt.rcParams.update({"font.family": "Arial", "font.size": 6, "axes.labelsize": 6.5,
                         "xtick.labelsize": 6, "ytick.labelsize": 5.5, "axes.titlesize": 7,
                         "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
                         "axes.linewidth": .5, "xtick.major.width": .5, "ytick.major.width": .5,
                         "xtick.major.size": 2.5, "ytick.major.size": 2.5,
                         "savefig.facecolor": "white", "figure.facecolor": "white"})
    # 183 x 200 mm, written as inches for source-level QA parsers.
    fig = plt.figure(figsize=(7.204724409448819, 7.874015748031496))
    fig.text(.04, .972, "CpG signature across adipocyte and twin adipose cohorts", fontsize=9, weight="bold")
    fig.text(.04, .943, "Unchanged source model  |  native hg19 gene context  |  independent computational replay", fontsize=6.5, color="#424950")
    ax = fig.add_axes([.245, .13, .305, .76])
    data = table[["discovery_delta_beta_pp", "replication_delta_beta_pp", "twin_delta_beta_pp"]].to_numpy()
    limit = float(np.ceil(np.nanmax(np.abs(data)) / 5) * 5)
    cmap = LinearSegmentedColormap.from_list("beta_direction", ["#336A8F", "#FFFFFF", "#B65758"])
    cmap.set_bad("#D5D9DD")
    mesh = ax.imshow(data, aspect="auto", cmap=cmap, norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit), interpolation="none")
    for row, column in np.argwhere(~np.isfinite(data)):
        ax.plot(column, row, marker="x", markersize=3, markeredgewidth=.5, color="#5C6268")
    labels = [f"{r.feature_id}  {r.gene_annotation_display}" for r in table.itertuples()]
    ax.set_yticks(np.arange(50), labels=labels)
    ax.tick_params(axis="y", length=0, pad=3)
    ax.set_xticks([0, 1, 2], labels=["Source", "EPIC", "Twins"])
    ax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False, length=0, pad=4)
    ax.set_xlim(-.5, 2.78)
    for index, hit in enumerate(table.twin_q_below_005_selected50):
        if hit:
            ax.plot(2.61, index, marker="D", markersize=2.7, color="#222222", clip_on=False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.text(.04, .918, "a", fontsize=8, weight="bold")
    fig.text(.245, .924, "Methylation difference", fontsize=7, weight="bold")
    cax = fig.add_axes([.275, .075, .235, .012])
    cb = fig.colorbar(mesh, cax=cax, orientation="horizontal", ticks=[-limit, 0, limit])
    cb.outline.set_linewidth(.5)
    cb.ax.tick_params(labelsize=6, length=2)
    cb.set_label("Delta beta (percentage points)", fontsize=6)
    fig.text(.045, .105, "50 CpGs; source coefficient rank", fontsize=5.5, color="#424950")
    fig.text(.60, .918, "b", fontsize=8, weight="bold")
    bx = fig.add_axes([.625, .60, .335, .285])
    bx.set_title("Discrimination and transfer", loc="left", pad=8, weight="bold")
    colors = ["#336A8F"] * 3 + ["#B65758"]
    for index, row in enumerate(aurocs.itertuples()):
        y = 3 - index
        bx.errorbar(row.auroc, y, xerr=[[row.auroc-row.ci_low], [row.ci_high-row.auroc]],
                    fmt="o", markersize=3.5, color=colors[index], elinewidth=.8, capsize=2, capthick=.7)
        bx.text(.47, y+.26, row.label, fontsize=5.8, va="bottom")
        bx.text(1.03, y-.22, f"{row.auroc:.3f} [{row.ci_low:.3f}, {row.ci_high:.3f}], n={row.n}",
                fontsize=5.6, va="top", ha="right", color="#424950")
    bx.axvline(.5, color="#B5BAC0", linewidth=.6, linestyle="--", zorder=0)
    bx.set(xlim=(.47, 1.03), ylim=(-.52, 3.6), xlabel="AUROC (conditional 95% CI)", yticks=[])
    bx.set_xticks([.5, .75, 1.0])
    bx.spines[["left", "right", "top"]].set_visible(False)
    fig.text(.625, .535, "Twin endpoint: higher versus lower BMI\nwithin discordant pairs; bulk SAT tissue.", fontsize=5.8, linespacing=1.4, color="#424950")
    fig.text(.60, .495, "c", fontsize=8, weight="bold")
    cx = fig.add_axes([.625, .24, .335, .22])
    cx.set_title("Paired twin loci, selected-50 FDR", loc="left", pad=8, weight="bold")
    hits = table.loc[table.twin_q_below_005_selected50].sort_values("q_value_selected_50_family", kind="stable")
    for index, row in enumerate(hits.itertuples()):
        y = len(hits)-index-1
        cx.errorbar(row.twin_delta_beta_pp, y,
                    xerr=[[row.twin_delta_beta_pp-row.twin_ci_low_pp], [row.twin_ci_high_pp-row.twin_delta_beta_pp]],
                    fmt="D", color="#B65758", markersize=3.5, elinewidth=.8, capsize=2, capthick=.7)
        region = " (body)" if row.gene_regions_native == "Body" else ""
        cx.text(-.3, y+.28, f"{row.feature_id}  {row.gene_annotation_display}{region}", fontsize=5.7, va="bottom")
        cx.text(6.55, y-.20, f"q={row.q_value_selected_50_family:.4f}", fontsize=5.6, ha="right", va="top", color="#424950")
    cx.axvline(0, linewidth=.6, linestyle="--", color="#B5BAC0", zorder=0)
    cx.set(xlim=(-.35, 6.55), ylim=(-.5, 3.6), xlabel="Higher minus lower BMI: delta beta (pp)", yticks=[])
    cx.set_xticks([0, 2, 4, 6])
    cx.spines[["left", "right", "top"]].set_visible(False)
    agree = int(table.twin_direction_matches_discovery.sum())
    unassigned = int(table.gene_annotation_display.eq("Unassigned").sum())
    fig.text(.625, .167, f"{agree}/50 directions agree; {len(hits)}/50 paired q < 0.05.\n{unassigned}/50 loci lack a native gene assignment.", fontsize=6.2, linespacing=1.45)
    fig.text(.625, .11, "Diamond: paired q < 0.05; grey x: unmeasured.\nPaired effects and CIs: 24 twin families.\nGene annotation does not imply expression causality.", fontsize=5.6, color="#424950", linespacing=1.4)
    # Fixed canvas retains the declared physical dimensions; no bbox expansion.
    fig.savefig(output / "obesity50_twin_gene_context.pdf", dpi=600)
    fig.savefig(output / "obesity50_twin_gene_context.svg", dpi=600)
    fig.savefig(output / "obesity50_twin_gene_context.png", dpi=600)
    fig.savefig(output / "obesity50_twin_gene_context.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)
    return {"directions_matching_discovery": agree, "paired_fdr_hits": len(hits), "unassigned_native_genes": unassigned,
            "EPIC_effect_unmeasured": int((~table.EPIC_effect_measured).sum()),
            "heatmap_color_limit_pp": limit, "statistical_family": "50 frozen-source-selected CpGs, not genome-wide"}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--portable-inputs", action="store_true", help="Read previously packaged inputs instead of repo receipts")
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    contract = output / "figure_contract.json"
    if not contract.exists():
        raise ValueError("A reviewed pre-render figure contract is required")
    inputs = output / "inputs"
    inputs.mkdir(exist_ok=True)
    sources = {"source_annotation_all50.tsv": ROOT / "docs/chronic_20261003/obesity_annotation_v1/cpg_gene_interpretation_all50.tsv",
               "twin_paired_gene_context.tsv": ROOT / "docs/twins_transfer_20261003_v2/paired_cpg_gene_context.tsv",
               "twin_summary.json": ROOT / "docs/twins_transfer_20261003_v2/summary.json",
               "twin_predictions.tsv": ROOT / "docs/twins_transfer_20261003_v2/all_predictions.tsv",
               "twin_source_audit.json": ROOT / "docs/twins_transfer_20261003_v2/source_audit.json",
               "obesity_replay_summary.json": ROOT / "docs/obesity_replay_20261003_v1/summary.json"}
    if not args.portable_inputs:
        for name, source in sources.items():
            shutil.copy2(source, inputs / name)
    table, aurocs = assemble(pd.read_csv(inputs / "source_annotation_all50.tsv", sep="\t"),
                            pd.read_csv(inputs / "twin_paired_gene_context.tsv", sep="\t"),
                            json.loads((inputs / "obesity_replay_summary.json").read_text()),
                            json.loads((inputs / "twin_summary.json").read_text()))
    table.to_csv(output / "annotated_all50_cpg.tsv", sep="\t", index=False)
    heatmap = table[["source_coefficient_rank", "feature_id", "gene_annotation_display", "gene_regions_native",
                     "discovery_delta_beta_pp", "replication_delta_beta_pp", "twin_delta_beta_pp", "twin_q_below_005_selected50"]]
    heatmap.to_csv(output / "panel_a_heatmap.tsv", sep="\t", index=False)
    aurocs.to_csv(output / "panel_b_discrimination.tsv", sep="\t", index=False)
    table.loc[table.twin_q_below_005_selected50].to_csv(output / "panel_c_paired_loci.tsv", sep="\t", index=False)
    checks = make_figure(table, aurocs, output)
    script_copy = output / "make_twins_cpg_figure.py"
    if Path(__file__).resolve() != script_copy.resolve():
        shutil.copy2(__file__, script_copy)
    receipt = {"generated_at": datetime.now(timezone.utc).isoformat(),
               "versions": {"python": platform.python_version(), "matplotlib": matplotlib.__version__, "numpy": np.__version__, "pandas": pd.__version__},
               "input_sha256": {name: sha(inputs / name) for name in sources}, "script_sha256": sha(script_copy),
               "contract_sha256": sha(contract), "scientific_data_checks": checks,
               "regenerate_portable": "python make_twins_cpg_figure.py --portable-inputs --output .",
               "limits": json.loads(contract.read_text())["limitations"]}
    (output / "regeneration_receipt.json").write_text(json.dumps(receipt, indent=2), encoding="utf-8")
    print(json.dumps(checks))


if __name__ == "__main__":
    main()

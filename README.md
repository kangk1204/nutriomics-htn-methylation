**English** | [한국어](README.ko.md)

# NutriOmics Disease Methylation Analysis

This repository classifies disease groups from DNA methylation data and links the selected CpG sites to genes and regulatory regions. The main analysis uses subcutaneous adipocytes from participants with and without obesity, with further evaluations in hypertension and twin adipose-tissue datasets.

[Full performance table](docs/phase1_extension_receipts_20261003/all_executed_extension_evaluations.tsv) | [Obesity analysis results](docs/obesity_replay_20261003_v1/summary.json) | [CpG and gene figure](docs/twins_figure_20261003_v1/obesity50_twin_gene_context.pdf) | [Project report, restricted access](https://github.com/kangk1204/nutriomics-final-report)

## Main results

These results come from analyses and reproduction runs completed in October 2026. AUROC measures how well model scores distinguish the two groups; values closer to 1 indicate better discrimination.

| Analysis | Dataset and participants | Evaluation | AUROC |
|---|---|---|---:|
| Obesity in subcutaneous adipocytes | GSE222595, 48 participants: 24 with obesity and 24 controls | Nested cross-validation, with feature selection and tuning repeated within the training folds | **0.932** |
| Separate obesity cohort | GSE67024, 29 women: 15 with obesity and 14 controls | Fixed 50-CpG model fitted on the training dataset | **0.981** |
| Obesity EPIC data | GSE222595, 48 participants: 24 with obesity and 24 controls | Cross-platform transfer of the model trained on 450K data | 0.915 |
| Twin adipose tissue | GSE68336, 24 monozygotic twin pairs, 48 participants | Discrimination between higher and lower BMI within families | 0.700 |
| Hypertension | GSE193795, 88 participants: 44 with hypertension and 44 controls | Elastic net result from the nested comparison of ridge, elastic net and linear SVM | 0.654 |

In the internal obesity evaluation, the number of features was selected within each training fold. The separate cohort received the fixed 50-CpG model fitted on the complete training dataset. The GSE67024 result reproduces the existing separate-cohort evaluation, with no retraining or tuning on that cohort. The EPIC analysis used the jointly preprocessed matrix from the original study.

The twin analysis evaluates within-family BMI differences in whole adipose tissue. Its tissue and outcome definition differ from the isolated-adipocyte obesity analysis. The [full performance table](docs/phase1_extension_receipts_20261003/all_executed_extension_evaluations.tsv) includes PR-AUC, Brier score, sensitivity, specificity and confidence intervals obtained by resampling the fixed predictions.

## CpG and gene interpretation

The 50 CpG sites are linked to reference coordinates, genes, transcripts and regulatory regions. Twin results include within-family methylation differences and multiple-testing corrections. Gene annotations describe the genomic location of each CpG; the direction of expression regulation is interpreted using separate supporting evidence.

![Gene annotations for obesity-associated CpGs and changes in twin adipose tissue](docs/twins_figure_20261003_v1/obesity50_twin_gene_context.png)

- [Annotations for all 50 CpGs](docs/twins_figure_20261003_v1/annotated_all50_cpg.tsv)
- [Within-family changes and gene context](docs/twins_transfer_20261003_v2/paired_cpg_gene_context.tsv)
- [Figure legend and analysis settings](docs/twins_figure_20261003_v1/legend.md)

## Analysis workflow

1. Read public methylation matrices and sample metadata, then confirm the disease and control groups.
2. Fit missing-value handling, feature selection, scaling and model tuning within each training fold.
3. Combine predictions for held-out samples to calculate internal performance.
4. Apply the fixed model to a separate cohort and save cohort-specific predictions and metrics.
5. Link CpG annotations and observed changes for biological interpretation.

The 44 participants with prehypertension are reserved for supplementary evaluation and excluded from the main hypertension training and tuning. Repeated measurements and twin analyses account for participant and family units. Each result file records the resampling unit and preprocessing settings.

## Results and reproduction

| Information | File |
|---|---|
| Obesity analysis settings | [protocol.json](docs/obesity_replay_20261003_v1/protocol.json) |
| Internal evaluation predictions | [source_repeated_oof.tsv](docs/obesity_replay_20261003_v1/source_repeated_oof.tsv) |
| Separate-cohort predictions | [gse67024_predictions.tsv](docs/obesity_replay_20261003_v1/gse67024_predictions.tsv) |
| Obesity reproduction script | [replay_obesity_primary.py](scripts/replay_obesity_primary.py) |
| Hypertension model comparison | [run_extension_htn_grid.py](scripts/run_extension_htn_grid.py) |
| Twin analysis | [evaluate_twins_transfer.py](scripts/evaluate_twins_transfer.py) |
| Initial hypertension evaluation and LOOCV | [Analysis summary](docs/validation/SUMMARY.md) |
| Training-set evaluation | [Separate diagnostic results](docs/trainset_20261003/SUMMARY_KO.md) |

Use Python 3.12 or later.

```sh
python -m pip install -e ".[test,annotation]"
nutriomics-methylation --help
python -m pytest -q
```

Run the basic hypertension analysis in the following order. Data downloads and training require disk space and computation time.

```sh
nutriomics-methylation fetch
nutriomics-methylation prepare
nutriomics-methylation train --cores 12 --spaces common full --methods nested5 loocv
nutriomics-methylation evaluate --bootstrap-replicates 1000
nutriomics-methylation export
```

For the obesity and twin analyses, use the scripts and analysis settings linked above. Source data and trained models are held in research data storage. This repository provides code, analysis settings, summary tables and predictions.

## Sources and licensing

The principal methylation datasets are NCBI GEO accessions GSE222595, GSE67024, GSE68336, GSE193795 and GSE42774. The repository name retains `htn` from the initial hypertension analysis; the current code also covers other disease evaluations.

The code is available under the [MIT License](LICENSE). External datasets retain their providers' terms, as described in [Software and data licenses](DATA_LICENSES.md).

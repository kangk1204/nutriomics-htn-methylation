# Chronic-disease extension: observed 3 October 2026 results

Seven official GEO accessions were acquired: GSE89632, GSE38642, GSE21232, GSE48325, GSE61258, GSE222595 and GSE135251. Six supplied usable numeric data, comprising seven source/target cohorts including the two separate obesity cohorts; GSE135251 was metadata-only in this analysis. Obesity uses the supplementary processed beta/detection-P matrix, since its series matrices contain no numeric rows. GSE38642 is expression, not methylation. Initial two truncated liver downloads were rejected; replacement files passed official byte length, SHA256 and full gzip CRC checks. Acquisition receipts retain both attempts.

All ten source cohort/ranking model results are retained: six initial screening combinations and four additional methylation source/target combinations. No outcome-driven model grid or target-cohort tuning is performed. Models use ridge logistic regression with balanced class weights, k10/50/200 and C0.1/1.0. F-statistic and variance ranking are separate prespecified alternatives. Missingness filtering, median imputation, feature selection and scaling are fitted inside every inner and outer training fold. Source nested5/minority-limited and inner nested3/minority-limited CV runs use seeds20261003/20261020/20261103. Repeated OOF AUROC is computed from the mean of three per-sample predictions, not the mean of three AUROCs.

| Cohort / modality | Controls / cases | F-statistic source OOF AUROC | Variance source OOF AUROC | Frozen target AUROC, F / variance |
|---|---:|---:|---:|---:|
| GSE89632, liver expression, NAFLD | 24 /39 | 0.966 | 0.960 | Not evaluated |
| GSE38642, islet expression, T2D | 54 /9 | 0.786 | 0.800 | Not evaluated |
| GSE21232, islet DNAm, T2D | 11 /5 | 0.800 | 0.945 | Not evaluated |
| GSE222595, SAT DNAm, obesity discovery | 24 /24 | 0.932 | 0.524 | 0.915 /0.627, replication24/24 |
| GSE48325, liver DNAm, NAFLD baseline | 34 /25 | 0.753 | 0.515 | 0.912 /0.412, GSE61258 controls32/cases21 |

[Initial screening](screen_v1/comparison.tsv) and [methylation transfer](methylation_transfer_v3/comparison.tsv) contain exact measured values. Their JSON summaries include AP, Brier, fixed0.5 sensitivity/specificity and confidence intervals. Prediction files, all150 outer-fold audits and inner split membership are versioned. [Verification](verification.json) independently recalculates48 prediction tables; [method tests](method_tests.xml) records19 passes. First failed parsing/serialization attempts and isolated test retry are described in [execution_attempts.json](execution_attempts.json); no methylation scores were produced in those failed attempts.

Obesity uses only SAT samples, one included SAT measurement per deposited sample accession, with native HA codes mapped through GEO description fields; the original study declares separate discovery/replication donor cohorts. Original-study190-sample QC is not reproduced; this analysis uses the deposited48SAT discovery and48SAT replication samples with detection-P masking. All192SAT/VAT assays are not counted as192independent people. The875,551 union probe rows are not all measured in each platform: discovery470,866 and replication844,023 have at least one finite value. [Measurement audit](prepared_measurement_audit.json) records this distinction and no entirely missing included sample. Author-processed data are used; upstream cohort normalization and raw-IDAT pipeline independence are not established.

GSE48325 postoperative and unknown disease-group samples are excluded before model fitting. GSE61258 excludes PSC/PBC and all DNAm-age-derived fields. Related recruitment and unknown donor linkage limit the two liver cohorts to exploratory transfer; absence of identical GSM/profile hashes does not prove independent donors. No individual age/sex adjustment, causal biomarker, dietary response, clinical performance or selected-positive confirmatory claim is made. Conditional bootstrap intervals use1,000 realized-prediction resamples and do not capture full model-selection/refit uncertainty. Screening several datasets remains exploratory.

## Reproduce

Use the versions in `screen_v1/receipt.json` and `methylation_transfer_v3/receipt.json`; acquire raw data into ignored `data/`, then verify its hashes against the preserved acquisition/source receipts before analysis. Use a fresh output directory for every real rerun. The names below describe the original complete run:

```bash
python scripts/acquire_chronic_geo.py --raw data/chronic_20261003/raw_v1 --docs docs/chronic_20261003/acquisition_v1
python scripts/download_chronic_assets.py --out data/chronic_20261003/assets_v1
python scripts/run_chronic_screen.py --raw data/chronic_20261003/raw_v1 --out docs/chronic_20261003/screen_v1 --cores 3
python scripts/run_chronic_methylation.py --assets data/chronic_20261003/assets_v1 --metadata docs/chronic_20261003/acquisition_v1 --prepared data/chronic_20261003/prepared_v3 --out docs/chronic_20261003/methylation_transfer_v3
python scripts/audit_chronic_prepared.py --prepared data/chronic_20261003/prepared_v3 --out docs/chronic_20261003/prepared_measurement_audit.json
python scripts/verify_chronic_outputs.py --repo . --out docs/chronic_20261003/verification.json
```

The source-acquisition receipts were copied from the CPU execution host. `collect_chronic_receipts.py` keeps native GEO metadata under ignored `data/` and renames versionable derived source receipts without changing bytes. To run new preparation, use newly acquired metadata at its original output location; do not expect the publication snapshot to contain raw arrays or copied native metadata. Prepared NPZ and trusted newly trained model artifacts remain only on the authorized CPU server and are excluded from Git.

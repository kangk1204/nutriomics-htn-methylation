# Claims supported by the new implementation

Verified public-data snapshot: 2026-10-01. This is separate from historical project data and historical slide claims.

| Quantity | Verified definition |
|---|---|
| 132 | Deposited GSE193795 GSM samples: HTN 44, control 44, preHT 44 |
| 88 | Development/OOF participants: HTN 44 + control 44 |
| 44 | preHT supplementary predictions, never fitted or tuned |
| 16 | GSE42774 exploratory external observations: 8 cases + 8 matched controls |
| 148 | Total public GSM inputs across the two GEO series, not newly recruited participants |
| 223,945 | CpGs in the deposited primary processed beta panel, not the complete 450k assay manifest |
| 26,486 | CpGs in the external processed beta panel |
| 10,749 | Measured CpG ID intersection frozen before phenotype mapping |
| 12 | Primary chip IDs, all containing HTN, preHT and control samples |
| 132 unique | Primary deposited descriptions; no identical beta profiles detected. Undisclosed relatedness remains unknown |

Observed Detection P failure fractions before fold-local filtering: primary 0.0161–0.1465%; external 0.0227–8.0722%. No sample was removed to improve an outcome metric. Fixed per-measurement Detection P masking is followed by training-fold missingness QC and imputation. External cohort averages/variances/outcomes never determine feature selection.

The estimand is discrimination of deposited hypertension cases from deposited controls. Case-control sampling does not identify population incidence, prospective risk or causal dietary response. User-facing outputs should be described as research model scores, not personal hypertension risk percentages.

## Historical HWP claims

| Original paragraph | Compatibility with this implementation |
|---|---|
| P0246, P0867–P0870: EPIC, 5 repeated P-value screens, Lasso, SHAP | These describe a different historical procedure. Newly implemented data are 450k/27k, training-fold variance ranking, ridge/elastic net, coefficient and sign stability. SHAP was not run. Do not describe this pipeline as an EPIC/SHAP revalidation |
| P0254, P0869, P0877: HTN 131/control 678, AUROC 0.88 | Source cohort, split, preprocessing and canonical performance evidence were not recovered. Do not reuse those denominators or AUROC as new results |
| P0253–P0254: 200 biomarkers, approximately 2,700 participants | Selected model variables and marker stability are exploratory. They are not newly validated clinical biomarkers and do not establish these historical totals |
| P0254/P0877: 2,000 cumulative samples and LOOCV milestone | New LOOCV applies to 88 development observations. The new public-data analysis cannot substantiate 2,000 historical recruited or cumulative specimens |
| P0247/P0254/P0264: December 2026 definitive totals | Current input evidence is dated October 1, 2026. Do not future-date the snapshot as a December measurement |
| P1579: platform superiority | No head-to-head clinical performance comparison with NutriGenomeDB, Eat4Genes or FoodAtlas was executed. Functional scope comparison is possible; superiority is unproven |

## Probe annotation audit

The official GPL13534 platform table within GSE193795 family SOFT confirms these reference annotations (genome coordinates as deposited, not newly remapped):

| CpG | Official gene field | CHR | MAPINFO | Region | In the primary processed beta panel? |
|---|---|---|---:|---|---|
| cg04983977 | GPR25 | 1 | 200842145 | TSS200 | Yes |
| cg15336997 | VCAN;VCAN;VCAN;VCAN | 5 | 82770250 | 5'UTR;5'UTR;5'UTR;5'UTR | No |
| cg20823859 | OVGP1 | 1 | 111971624 | TSS1500 | No |

Annotation agreement does not establish historical feature importance, replicate the Circulation OVGP1 result, or validate a clinical marker. The processed panel excludes the last two loci even though the complete deposited platform annotation contains them. The reason for the original processed-panel filtering is not fully described in the brief deposited processing metadata; raw IDAT reprocessing would be a separate analysis.

## Validation limits

All candidate models and both feature spaces are reported. The best external result must not be selected as an unbiased primary estimate. AUROC/PR/Brier intervals condition on realized OOF predictions; they do not repeat the entire nested fitting experiment. Nested CV still shares overlapping training sets. Clinical biomarkers and clinical calibration require larger independent cohorts with individual covariates and outcomes.

Random participant folds share some technical chips, although the chip-by-phenotype table is mixed rather than completely confounded. New-chip performance is not estimated by those internal folds. No blanket ComBat, joint external normalization or invented age/sex adjustment was applied. The small external cohort probes transfer to a different assay/population but cannot isolate the source of domain shift.

DNA methylation is not gene-expression log fold change. Promoter proximity does not guarantee an inverse expression relationship. Further integration must preserve CpG-to-gene ambiguity, require an independently supported eQTM or treat pathway relations as unvalidated mechanistic hypotheses.

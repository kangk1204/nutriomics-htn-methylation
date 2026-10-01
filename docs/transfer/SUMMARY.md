# Exploratory reciprocal transfer

Supplementary direction-of-transfer experiment. Previously evaluated target outcomes preclude a new untouched external validation claim; no improvement claim from selecting the better direction.

Source: GSE42774 n16 (8 published matched pairs). Target: previously evaluated GSE193795 HTN44/control44; preHT44 excluded.
All preprocessing and hyperparameters are learned from source training groups only. No target normalization, batch correction, feature selection, threshold selection or coefficient flipping.

| Analysis | N | AUROC (95% conditional interval) | Brier | Sensitivity | Specificity |
|---|---:|---|---:|---:|---:|
| ridge/source_nested4_pair_cv | 16 | 0.641 (0.297–0.922) | 0.233 | 0.750 | 0.625 |
| ridge/target_reciprocal_transfer | 88 | 0.570 (0.451–0.693) | 0.270 | 0.841 | 0.273 |
| elasticnet/source_nested4_pair_cv | 16 | 0.516 (0.187–0.844) | 0.284 | 0.500 | 0.500 |
| elasticnet/target_reciprocal_transfer | 88 | 0.562 (0.448–0.679) | 0.276 | 0.568 | 0.409 |

Intervals condition on fixed realized predictions. Tiny source size and cross-platform/population differences remain. No individual age/sex adjustment was invented.
Original forward-transfer and primary nested-CV results remain frozen and must be reported beside these supplementary results. A more favorable point estimate does not establish improved generalization.

# Observed validation results

Research-only results on deposited public samples. No guaranteed performance, clinical risk or causal dietary response.

All intervals condition on the realized predictions. They do not include CV/model-selection uncertainty.

| Analysis | N | AUROC (95% conditional interval) | Average precision | Brier |
|---|---:|---|---:|---:|
| common_elasticnet/external | 16 | 0.641 (0.312–0.922) | 0.611 | 0.254 |
| common_elasticnet/loocv | 88 | 0.470 (0.353–0.600) | 0.492 | 0.257 |
| common_elasticnet/nested5 | 88 | 0.541 (0.426–0.650) | 0.539 | 0.291 |
| common_ridge/external | 16 | 0.562 (0.219–0.875) | 0.611 | 0.257 |
| common_ridge/loocv | 88 | 0.468 (0.346–0.590) | 0.507 | 0.286 |
| common_ridge/nested5 | 88 | 0.468 (0.343–0.580) | 0.481 | 0.305 |
| full_elasticnet/loocv | 88 | 0.554 (0.433–0.684) | 0.540 | 0.268 |
| full_elasticnet/nested5 | 88 | 0.569 (0.448–0.696) | 0.558 | 0.283 |
| full_ridge/loocv | 88 | 0.366 (0.250–0.489) | 0.436 | 0.392 |
| full_ridge/nested5 | 88 | 0.432 (0.312–0.560) | 0.454 | 0.303 |

The common probe mask was frozen from measurement IDs before outcome mapping. All learned QC, median imputation, variance selection, scaling and tuning used training folds only.

preHT44 was excluded from fit/tuning. GSE42774 is an exploratory n16 transfer study across 450k/27k, age, ancestry and sex distributions. Individual age/sex adjustment was not fabricated.

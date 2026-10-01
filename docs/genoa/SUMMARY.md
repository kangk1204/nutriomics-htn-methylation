# Independently deposited GENOA measured-BP proxy evaluation

Source88 common-probe ridge and elastic-net models remain frozen. Target endpoint is measured-BP proxy, not clinical essential hypertension.
No target-cohort fitting/tuning, covariate adjustment, threshold selection or prediction inversion. Unknown family dependence and medication correction remain.
Coverage-limited deployment: 10689/10749 source probe IDs measured; 60 unprovided probes stay NaN and use source-fitted medians. No complete platform-coverage claim.

| Analysis | N | AUROC (95% conditional interval) | Brier | Sensitivity | Specificity |
|---|---:|---|---:|---:|---:|
| ridge/all_proxy | 179 | 0.551 (0.453–0.641) | 0.252 | 0.516 | 0.604 |
| ridge/age_50_65_secondary | 67 | 0.488 (0.343–0.633) | 0.267 | 0.381 | 0.560 |
| elasticnet/all_proxy | 179 | 0.516 (0.410–0.614) | 0.234 | 0.968 | 0.057 |
| elasticnet/age_50_65_secondary | 67 | 0.545 (0.398–0.699) | 0.241 | 1.000 | 0.120 |

Case SBP≥140 or DBP≥90; control SBP<120 and DBP<80; middle BP excluded. Age50–65 is a prespecified secondary descriptive population, not matching/covariate adjustment.
Participant-bootstrap intervals condition on fixed predictions and do not account for unknown family dependence or model-selection variability. All original clinical-HTN and reciprocal-transfer results must remain reported separately.

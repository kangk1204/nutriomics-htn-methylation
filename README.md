# NutriOmics hypertension methylation

A research implementation using real public peripheral-blood DNA methylation data. This package separates a disease classifier from dietary-response evidence. Its probabilities are not clinical risk estimates, and methylation coefficients do not establish expression direction or causal dietary effects.

The actual frozen validation run is complete. [Observed results](docs/validation/SUMMARY.md) and the [Korean conclusion](docs/validation/CONCLUSION_KO.md) show that the specified models did not establish useful hypertension discrimination. The main nested 5-fold elastic-net AUROCs were 0.541 (common panel) and 0.569 (full processed panel). External models classified all 16 samples as positive at the fixed 0.5 threshold. Historical AUROC 0.88 is not a result of this implementation.

The separate [3 October full-training diagnostic](docs/trainset_20261003/SUMMARY_KO.md) refits the same88 development samples and scores those exact samples: AUROC0.902–0.998. These are training resubstitution scores; original CV/transfer results are preserved. The user's small-sample exploratory development context is recorded without reconstructing an unavailable historical model.

The [chronic-disease extension](docs/chronic_20261003/README.md) acquires seven GEO accessions and evaluates all ten declared cohort/ranking combinations. NAFLD liver expression gives repeated nested OOF AUROC0.966; obesity subcutaneous-adipocyte DNAm gives source OOF0.932 and frozen450K→EPIC cohort AUROC0.915 (48source/48target samples). The author-processed combined matrix's upstream normalization was not independently reconstructed. Lower variance-only results and exploratory liver transfer are retained. These tissue classifiers do not establish blood-based clinical or dietary effects.

## Data contract

* [GSE193795](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193795): processed 450k measurements, HTN 44/control 44 for development. preHT 44 is excluded from every fit, feature selection and hyperparameter search and used only for supplementary predictions.
* [GSE42774](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE42774): independent 27k, HTN 8/control 8, young African American males aged 14–23, with eight deposited matched pairs. This is exploratory transfer across platform and population, not established clinical generalization.
* Official GEO HTTPS series matrices supply beta values and sample metadata. Official series-family SOFT files supply observed Detection P values and platform gene annotations. The package does not need the enormous platform-family files or the 1.2 GB primary IDAT archive.
* The downloaded primary processed panel contains 223,945 CpGs; it is not all probes on the 450k array. Filtering performed by the original depositors cannot be undone here. The measured primary/external intersection contains 10,749 CpGs, frozen from IDs before mapping outcomes. Full-panel models and common-panel transferable models are reported separately.
* Detection P > 0.01 and nonfinite Detection P mask individual measurements using the same prespecified rule in each dataset. No external cohort averages, variances or outcomes determine feature QC. Missingness filtering, median imputation, variance ranking, scaling and tuning are learned within training folds.
* Deposited primary individual age and sex are unavailable. External sex is deposited as Male; individual age is unavailable. Age/sex adjustment is not fabricated from publication-level demographics. Twelve primary chips each contain all three phenotypes; the chip-by-phenotype table is exported. Unique deposited accessions and duplicate beta profiles are checked; undisclosed familial relationships cannot be ruled out.

## Reproduce

Python 3.12 or later. CPU only; at most 12 worker processes, one numerical-library thread per worker. No GPU or H200 is used.

```bash
python -m pip install -e '.[test]'
python -m pytest -q
nutriomics-methylation fetch
nutriomics-methylation prepare
nutriomics-methylation train --cores 12 --spaces common full --methods nested5 loocv
nutriomics-methylation evaluate --bootstrap-replicates 1000
nutriomics-methylation export
```

The CLI also accepts `--root /absolute/project/path`. Run `scripts/run_remote.sh` on the designated CPU server after transferring the prepared input bundle.

Primary evaluation is nested stratified 5-fold, inner stratified 5-fold. Nested outer LOOCV is a supplementary comparison. Every inner fold fits its own missingness filter (training missing fraction ≤ 0.10), medians, variance ranking and scaler. The inner grid is k={100,500}, C={0.01,0.1,1}; elastic-net l1 ratio={0.1,0.5}. Ridge uses liblinear and elastic net uses saga. Selection maximizes mean inner AUROC; ties retain the first grid candidate. A final inner-tuned refit uses development samples only. External outcomes never choose a model, threshold or features. Both model families are reported regardless of external performance.

The threshold is prespecified 0.5. Outputs include pooled OOF AUROC, average precision, trapezoidal PR area, Brier, calibration slope/intercept, sensitivity/specificity, confidence intervals, fold audits and marker stability. AUROC is never computed on a single LOOCV test observation. Discrimination/Brier intervals bootstrap realized participant predictions (matched pairs for external); they do not refit models or capture model-selection uncertainty. Sensitivity/specificity intervals are exact binomial. preHT has no fabricated binary truth; only its probability distribution is summarized.

## Outputs and provenance

An additional reciprocal experiment is available as `nutriomics-methylation transfer --root /absolute/project/path --cores 12`. It trains on GSE42774 n16 and scores the previously evaluated GSE193795 HTN44/control44, with matched-pair source CV and source-only preprocessing/tuning. It writes separate `results_transfer/` and `docs/transfer/` outputs. [Protocol](docs/RECIPROCAL_TRANSFER_PROTOCOL.md) and [observed supplementary results](docs/transfer/SUMMARY.md) explain why this is exploratory and cannot establish improvement by choosing the more favorable transfer direction. Observed target AUROC was ridge0.570 / elastic net0.562; both conditional intervals include0.5. Frozen primary results remain unchanged.

`data/source_manifest.json` records official URLs, sizes, SHA256 and acquisition checks. `data/prepared/` contains label-aligned NPZ, sample metadata, detection QC, fixed probe mask, chip contingency table and official gene annotations. Raw inputs are excluded from Git. `results/` contains fold membership, inner candidate scores/convergence warnings, OOF predictions, fitted pipelines, untouched preHT predictions, exploratory external predictions, all-marker stability tables and SVG/PNG plots. `docs/validation/` is a small versionable export of real results and source/artifact hashes. No synthetic benchmark is presented as a real result.

Changing input bytes, grid, seed, feature spaces or runtime versions while resuming an existing results directory fails rather than silently mixing runs. Clear/move a result directory deliberately before beginning a different analysis. A complete marker export is exploratory; correlated CpGs and variance filtering can yield unstable genes. Gene annotations do not resolve gene-expression or target effect direction.

## Primary scientific sources

* Bai etal. Circulation2022. [DOI10.1161/CIRCULATIONAHA.121.057178](https://doi.org/10.1161/CIRCULATIONAHA.121.057178), [PMID36172862](https://pubmed.ncbi.nlm.nih.gov/36172862/): primary study associated with GSE193795.
* AlfonsoPerez and DelgadoMartinez. J.PersonalizedMedicine2023. [DOI10.3390/jpm13050787](https://doi.org/10.3390/jpm13050787), [PMID37240957](https://pubmed.ncbi.nlm.nih.gov/37240957/): an existing predictive reanalysis of the same primary cohort; its accuracy is not a result of this package.
* Wang etal. PLoSONE2013. [DOI10.1371/journal.pone.0053938](https://doi.org/10.1371/journal.pone.0053938), [PMID23325143](https://pubmed.ncbi.nlm.nih.gov/23325143/): external genome-wide n16 cohort. Additional targeted-pyrosequencing participants are not genome-wide test samples.

MIT applies to newly implemented code. Source datasets and third-party resources retain their own terms. No external code has been copied into this implementation.

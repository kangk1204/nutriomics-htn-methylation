# Frozen obesity classifier CpG interpretation, 2026-10-03

The interpretation uses the actual 50 CpGs and coefficients of the source-only
frozen classifier, not a new model fit. Official GPL13534 self-only SOFT was fully
read (485,577 rows, gzip CRC pass), and all selected IDs mapped exactly once.
Native coordinates are hg19/build37, native MAPINFO rather than zero-based BED.

39 of 50 probes have native RefGene names, representing 22 unique names. Transcript
duplicates are collapsed only for display; full accessions, context pairs and
native symbols are retained. Eleven probes have no native gene assignment.
Native names are versioned annotations, not necessarily current HGNC symbols.

`cpg_gene_interpretation_all50.tsv` contains all exact probes, standardized source
coefficients, observed beta differences and conditional descriptive intervals,
selection frequencies across 15 source outer fits, native SNP flags and distinct
published annotation layers. Top 12 are chosen by absolute source coefficient.
All 46 probes observed in both cohorts have the same mean-difference direction.
Four target-unobserved probes remain NA, not assigned an invented replication
effect; the already frozen prediction uses source-learned median imputation.

`published_exact_probe_matches.json` is a reformatted subset of McAllan et al.,
Nature Communications 14, 2784 (2023), DOI
[10.1038/s41467-023-38439-z](https://www.nature.com/articles/s41467-023-38439-z).
Exact CpG joins find 22 ST2 sentinels, 3 probes in ST7 (4 genic expression links),
1 in ST9 (2 functional-mapping/expression links), and 4 in ST11 (7 shared-TAD
expression links). Original RNA-seq was not downloaded/reanalysed. These are
reused same-study associations, not independent biological replication or
experimental causal validation.

Original article/supplement: CC-BY-4.0
([license](https://creativecommons.org/licenses/by/4.0/)); attribution preserved,
selected rows reformatted, no original workbook redistributed. Official platform:
[GPL13534](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL13534).
Acquisition receipts record URL, source hash, size and full-read status. Raw
downloads and participant matrices remain in ignored data/, model in models/.

Reproduce on the CPU host holding the already prepared matrices and frozen model:

```sh
python -m pip install -e '.[annotation]'
python scripts/acquire_obesity_annotation.py .
python scripts/interpret_obesity_cpgs.py .
```

Do not use these procedures to replace the frozen performance results. Group
differences are post hoc, unadjusted for individual age/sex and not selection-
adjusted EWAS tests. Eight native SNP/probe flags and cross-reactivity review
limitations are disclosed. Genomic overlap, nearest gene, functional mapping,
shared TAD, model weight and expression association are separate quantities.

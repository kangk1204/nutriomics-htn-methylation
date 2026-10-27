# Software and data licenses

## Project software

The original code and documentation in this repository are distributed under the [MIT License](LICENSE). Dependencies installed through Python, R or Bioconductor retain their own licenses.

## Research data

The software license applies to this project's implementation. External datasets, annotations, publications and pretrained model weights retain the terms set by their providers. Cite the original study and accession when using a research dataset. Input hashes, analysis settings and source identifiers accompany the results in `docs/`.

| Resource | Use in this project | Source and terms |
|---|---|---|
| NCBI GEO | Public methylation and expression studies, sample annotations and processed subsets | [GEO disclaimer](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html) and [GEO FAQ](https://www.ncbi.nlm.nih.gov/geo/info/faq.html). NCBI places no restrictions on use or distribution of GEO data; rights held by submitters still apply. Cite the accession and original publication. |
| BioStudies S-BSST831 | Endocrine hypertension analysis | [Study record](https://www.ebi.ac.uk/biostudies/studies/S-BSST831) and [original publication](https://pmc.ncbi.nlm.nih.gov/articles/PMC9635165/). Study and sample identifiers are retained in the analysis records. |
| GEO platform annotations | CpG coordinates, gene context and probe flags | Platform accessions, including GPL13534, identify the annotation source. Manufacturer and submitter notices remain applicable to these annotations. |

The principal obesity analyses use GSE222595, GSE67024 and GSE68336. The initial hypertension analyses use GSE193795 and GSE42774. Additional studies are identified in their corresponding result directories.

Files in `docs/` include predictions, selected CpG measurements, sample mappings, annotations and figures derived from these public studies. Sample identifiers refer to the original public records. Figure source tables and analysis scripts are supplied alongside the results. The MIT software license does not change the rights in the underlying study data or platform annotations.

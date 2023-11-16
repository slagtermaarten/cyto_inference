# Cytokine sensing inference for single cells using (single cell) transcriptome data

This repo includes all of the code used in the publication by Hoekstra &amp; Slagter et al.: 'Distinct spatiotemporal dynamics of CD8+ T cell-derived cytokines in the tumor microenvironment'.
This project involved the use of single-cell transcriptomic read-outs to infer the dissemination of T cell released cytokines through the tumor microenvironment.

This is a view into the code that preprocessed the raw data and generated the figures for use in the manuscript. 
The raw data is available on [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220738). 
The associated preprocessed data is available on [Mendeley data](https://doi.org/10.17632/2wwjdppm7f.2).

All analyses were done in `R` 4.1. Preprocessing is done using the `targets` pipelining system, rules are included in the `targets` directory.

This repo assumes its root directory to be sym-linked to `~/MirjamHoekstra`. 
A set of R variables is asssumed to be encoded in an `.Renviron` file located under the project roor, of which I include here a modified/mock version for your reference.


# cyto_inference

This repo includes all of the code used in a manuscript by Hoekstra &amp; Slagter et al., which involves the use of single-cell transcriptomic read-outs to infer the spreading of cytokines through the tumor microenvironment.

This is currently *merely* a view into the code that preprocessed the raw data and generated the figures for use in the manuscript. As the raw source data isn't publically accessible (yet), others won't be able to replicate the (pre-)processed objects relying on these raw data, nor the figures listed in the Rmarkdown files under `rmd`. This might change in the future.

All analyses were done in `R` 4.1. Preprocessing is done using the `targets` pipelining system, rules are included in the `targets` directory.

This repo assumes its root directory to be sym-linked to `~/MirjamHoekstra`. A set of R variables is asssumed to be encoded in an .Renviron file, that I did not include in this version of the repo for security reasons.

# Cytokine sensing inference for single cells using (single cell) transcriptome data

This repo includes all of the code used in the publication by Hoekstra &amp; Slagter et al.: 'Distinct spatiotemporal dynamics of CD8+ T cell-derived cytokines in the tumor microenvironment'.
This project involved the use of single-cell transcriptomic read-outs to infer the dissemination of T cell released cytokines through the tumor microenvironment.

This is a view into the code that preprocessed the raw data and generated the figures for use in the manuscript. 
The raw data is available on [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE220738). 
The associated preprocessed data is available on [Mendeley data](https://doi.org/10.17632/2wwjdppm7f.2).

All analyses were done in `R` 4.1. Preprocessing is done using the `targets` pipelining system, rules are included in the `targets` directory.

This repo assumes its root directory to be sym-linked to `~/MirjamHoekstra`. 
A set of R variables is asssumed to be encoded in an `.Renviron` file located under the project root, of which I include here a modified/mock version for your reference.

General reproduction instructions (untested)
* Install a version R.4.1 or up on a Unix-like system
* git clone this repo, e.g. `git clone https://github.com/slagtermaarten/cyto_inference.git`, on the system from step 1
* (Sym)link the resulting directory to `~/MirjamHoekstra`: `ln -s ./cyto_inference ~/MirjamHoekstra`
* Change directory to `cyto_inference`/`~/MirjamHoekstra`, install the project's dependencies by running `Rscript bin/install_deps.R`. `renv` and I have never gotten to be friends.
* Update the `.Renviron` configuration file to reflect your preferences and file system
* Download the preprocessed (single cell) RNASeq data from Mendeley and place the resulting directories in a dir called `data`. Rename the experiment dirs and add the prefix `raw_` to each of them, i.e. `cd data; find ./ -name 'exp' -type d -maxdepth 1 -exec mv {} raw_{} \;`
* Go back to the project root and run `Rscript bin/tar_make` to populate the `targets` object store
* You should now be ready to explore the notebooks and access the `targets` objects that are referenced within them

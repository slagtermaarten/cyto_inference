options(Ncpus = 32)
R_lib_dir <- Sys.getenv("R_LIBS_USER")
.libPaths(R_lib_dir)

install_github <- function(s) {
  if (!requireNamespace('remotes', quietly = TRUE))
      install.packages('remotes')
  remotes::install_github(s)
}

install <- function(s) {
  if (!requireNamespace('devtools', quietly = TRUE))
      install.packages('devtools')
  devtools::install(s)
}

try_install <- function(p) {
  R_lib_dir <- Sys.getenv("R_LIBS_USER")
  if (!dir.exists(R_lib_dir)) dir.create(R_lib_dir)
  if (!requireNamespace('BiocManager', quietly = TRUE))
      install.packages('BiocManager')
  if (!require(p, character.only = T, quietly = T)) {
    BiocManager::install(p, ask = F, lib = R_lib_dir)
  }
}


cit_pkgs <- c('devtools', 'ComplexHeatmap', 'credentials', 'milo', 'sctransform', 'voom', 'glmnet')
pkgs <- c('miloR', 'car',
  'parsnip', 'tidymodels', 'yardstick', 'ranger', 'kknn', 
  'devtools', 'curl', 'targets', 'ggpubr', 'ArrayExpress',
  'umap', 'progressr', 'rsample', 
  'Hmisc', 
  'glmnet',
  'preprocessCore', 'envnames', 'voom',
  'pracma', 'targets', 
  'DESeq2', 'wesanderson',
  'biomaRt', 'furrr', 'fgsea', 'naturalsort',
  'tidyverse', 'knitr', 
  'tikzDevice',
  'scran', 'SingleCellExperiment', 'limma', 'edgeR', 'ROCR',
  'devtools', 'Seurat', 'GGally', 'knitr', 'pacman', 'circlize',
  'NMF', 'e1071')
already_cited <- c('targets', 'glmnet', 'milo', 'sctransform', 'voom')


for (p in pkgs) {
  try_install(p)
}

# U_lib_dir <- if (dir.exists('~/libs')) '~/libs' else '/u_libs'
if (F && !require('credentials', character.only = T, quietly = T)) {
  remotes::install_github("r-lib/credentials", lib = R_lib_dir)
}

if (!require('ComplexHeatmap', character.only = T, quietly = T)) {
  install_github('Jokergoo/ComplexHeatmap')
  library(ComplexHeatmap)
}

if (!require('maartenutils', character.only = T, quietly = T)) {
  install_github('slagtermaarten/maartenutils')
}

if (!require('genesets', character.only = T, quietly = T)) {
  install_github('slagtermaarten/genesets')
}

source('~/MirjamHoekstra/R/init.R')

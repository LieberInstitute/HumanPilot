library('SingleCellExperiment')
library('here')
library('jaffelab')
library('scater')
library('scran')
library('pheatmap')
library('readxl')
library('Polychrome')
library('cluster')
library('limma')
library('sessioninfo')
library('limma')

load("rda/sce_layer.Rdata", verbose = TRUE)
source("layer_specificity_functions.R", echo = TRUE)

## Use genes_km_raw from layer_specificity_functions.R
## Match the layer names, combine 6 and 6b into one

genes_km_mat <- genes_km_raw[, c(1:6, 'WM')]
colnames(genes_km_mat)[-7] <- paste0('Layer', colnames(genes_km_mat)[-7])
genes_km_mat[['Layer6']] <- pmin(genes_km_raw[['6']] + genes_km_raw[['6b']], 1)



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

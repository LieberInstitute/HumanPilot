###
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
library('qlcMatrix')

## Load data
load(here(
    'Analysis',
    'Human_DLPFC_Visium_processedData_sce_scran.Rdata'
))

# get expression
mat = assays(sce)$logcounts

## filter
exprsIndex = rowMeans(mat) > 0
mat = mat[exprsIndex,]

## calc dist
dd = dist(t(mat))


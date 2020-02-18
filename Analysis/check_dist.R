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

## Load data
load(here(
    'Analysis',
    'Human_DLPFC_Visium_processedData_sce_scran.Rdata'
))

## filter to variable genes
sce_hvg = sce[top.hvgs,]

## calc dist
dd = dist(t(assays(sce_hvg)$logcounts))


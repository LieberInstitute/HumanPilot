##
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
library('reshape2')
library('lmerTest')
library('WGCNA')

## multithread
allowWGCNAThreads(6)

## Load data
load(here(
    'Analysis',
    'Human_DLPFC_Visium_processedData_sce_scran.Rdata'
))

## load layer guesses
load(here("Analysis", "Layer_Guesses", "rda", "layer_guess_tab.Rdata"))

## add layer info
layer_guess_tab$layer[layer_guess_tab$layer == 'Layer 2/3'] <- 'Layer 3'
layer_guess_tab$Layer = gsub("ayer ", "", layer_guess_tab$layer)

## add
sce$Layer = layer_guess_tab$Layer[match(sce$key, layer_guess_tab$key)]

## remove all 0s
gIndex = rowSums(assays(sce)$logcounts) > 0
sce = sce[gIndex,]

#####################
## split by sample ##
#####################

sIndexes = splitit(sce$sample_name)
sce_list = lapply(sIndexes, function(ii) sce[,ii])
exprs_list = lapply(sce_list, function(x) assays(x)$logcounts)

###################
## WGCNA power ####
###################

powers <- c(1:10, seq(from = 12, to=20, by=2))
sftthresh_list <- mclapply(exprs_list, function(geneExprs) {
		pickSoftThreshold(t(geneExprs), powerVector = powers,
           networkType = "signed", verbose = 5)
}, mc.cores=6)
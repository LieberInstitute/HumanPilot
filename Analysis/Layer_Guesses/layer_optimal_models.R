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

## Use genes_km_raw and genes_bm_raw from layer_specificity_functions.R
## Match the layer names, combine 6 and 6b into a single layer for thje KM list
## Combine 5a and 6b for the BM list
genes_km_mat <- genes_km_raw[, c(1:6, 'WM')]
colnames(genes_km_mat)[-7] <- paste0('Layer', colnames(genes_km_mat)[-7])
genes_km_mat[['Layer6']] <- pmin(genes_km_raw[['6']] + genes_km_raw[['6b']], 1)

genes_bm_mat <- genes_bm_raw[, as.character(c(1:4))]
genes_bm_mat[['5']] <- pmin(genes_bm_raw[['5a']] + genes_bm_raw[['5b']], 1)
genes_bm_mat[['6']] <- genes_bm_raw[['6']]
colnames(genes_bm_mat) <- paste0('Layer', colnames(genes_bm_mat))
genes_bm_mat[['WM']] <- 0

## Genes to query
genes_query <- data.frame(
    symbol = c(genes_km_raw$Gene, genes_bm_raw$Gene),
    list = rep(c('KM', 'BM'), c(nrow(genes_km_raw), nrow(genes_bm_raw))),
    stringsAsFactors = FALSE
)
genes_query$m <- match(tolower(genes_query$symbol), tolower(rowData(sce_layer)$gene_name))
genes_query$gene_name <- rowData(sce_layer)$gene_name[genes_query$m]
genes_query$ensembl <- rownames(sce_layer)[genes_query$m]

## Combine the matrices into one
genes_query_mat <- rbind(
    genes_km_mat, genes_bm_mat
)

## Keep only those that are present in the data and have useful matrix info
## otherwise the model fails later
marker_genes <- list(
    gene_info = genes_query[!is.na(genes_query$m) & rowSums(genes_query_mat) > 0, ],
    layer_mat = genes_query_mat[!is.na(genes_query$m) & rowSums(genes_query_mat) > 0, ]
)



##### Compute the duplicate correlation globally
## Extract the data
mat <- assays(sce_layer)$logcounts

## Build a group model
mod <- with(colData(sce_layer), model.matrix( ~ 0 + layer_guess))
colnames(mod) <- gsub('layer_guess', '', colnames(mod))
## Takes like 2 min to run
corfit <-
    duplicateCorrelation(mat, mod, block = sce_layer$subject_position)
#######


## Process each marker gene
layer_idx <- splitit(sce_layer$layer_guess)
eb_markers_list <- lapply(seq_len(nrow(marker_genes$gene_info)), function(i) {
    res <- rep(0, ncol(sce_layer))
    layers <- colnames(marker_genes$layer_mat)[marker_genes$layer_mat[i, ] > 0]
    res[unlist(layer_idx[layers])] <- 1
    m <- model.matrix( ~ res)
    eBayes(
        lmFit(
            mat[marker_genes$gene_info$m[i], , drop = FALSE],
            design = m,
            block = sce_layer$subject_position,
            correlation = corfit$consensus.correlation
        )
    )
})


pvals_markers <- sapply(eb_markers_list, function(x) {
    x$p.value[, 2, drop = FALSE]
})
colSums(data.frame(
    'FDRsig' = p.adjust(pvals_markers, 'fdr') < 0.05,
    'Pval10-6sig' = pvals_markers < 1e-6,
    'Pval10-8sig' = pvals_markers < 1e-8
))
# FDRsig Pval10.6sig Pval10.8sig 
#     87          50          36 
length(pvals_markers)
# [1] 130

tstats_markers <- sapply(eb_markers_list, function(x) {
    x$t[, 2, drop = FALSE]
})

## Create a summary
marker_genes_summary <- cbind(marker_genes[[1]], marker_genes[[2]])
marker_genes_summary$p_value <- pvals_markers
marker_genes_summary$fdr <- p.adjust(pvals_markers, 'fdr')
marker_genes_summary$t_stat <- tstats_markers


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

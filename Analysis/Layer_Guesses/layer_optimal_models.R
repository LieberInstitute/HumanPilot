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
colnames(genes_km_mat)[-7] <-
    paste0('Layer', colnames(genes_km_mat)[-7])
genes_km_mat[['Layer6']] <-
    pmin(genes_km_raw[['6']] + genes_km_raw[['6b']], 1)

genes_bm_mat <- genes_bm_raw[, as.character(c(1:4))]
genes_bm_mat[['5']] <-
    pmin(genes_bm_raw[['5a']] + genes_bm_raw[['5b']], 1)
genes_bm_mat[['6']] <- genes_bm_raw[['6']]
colnames(genes_bm_mat) <- paste0('Layer', colnames(genes_bm_mat))
genes_bm_mat[['WM']] <- 0

## Genes to query
genes_query <- data.frame(
    symbol = c(genes_km_raw$Gene, genes_bm_raw$Gene),
    list = rep(c('KM', 'BM'), c(nrow(genes_km_raw), nrow(genes_bm_raw))),
    stringsAsFactors = FALSE
)
genes_query$m <-
    match(tolower(genes_query$symbol), tolower(rowData(sce_layer)$gene_name))
genes_query$gene_name <- rowData(sce_layer)$gene_name[genes_query$m]
genes_query$ensembl <- rownames(sce_layer)[genes_query$m]

## Combine the matrices into one
genes_query_mat <- rbind(genes_km_mat, genes_bm_mat)

## Keep only those that are present in the data and have useful matrix info
## otherwise the model fails later
marker_genes <- list(gene_info = genes_query[!is.na(genes_query$m) &
        rowSums(genes_query_mat) > 0, ],
    layer_mat = genes_query_mat[!is.na(genes_query$m) &
            rowSums(genes_query_mat) > 0, ])

## How many of these are unique models?
marker_genes$gene_info$models <-
    apply(marker_genes$layer_mat > 0, 1, function(x) {
        paste(colnames(marker_genes$layer_mat[x]), collapse = '+')
    })
length(unique(marker_genes$gene_info$models))
# [1] 29

## Check if any genes are duplicated and if the layer assignment agree
dups_ensembl <-
    unique(marker_genes$gene_info$ensembl[duplicated(marker_genes$gene_info$ensembl)])
dups_models <- lapply(dups_ensembl, function(ens) {
    i <- which(marker_genes$gene_info$ensembl == ens)
    marker_genes$gene_info$models[i]
})
names(dups_models) <- dups_ensembl

## At most a gene shows up twice: once per list
stopifnot(all(lengths(dups_models) == 2))

## Find which ones are truly duplicated: the gene and the model are the same
dups_ensembl_remove <-
    dups_ensembl[lengths(sapply(dups_models, unique)) == 1]

## Update the list section to label them as fully duplicated in KM and BM
marker_genes$gene_info$list[marker_genes$gene_info$ensembl %in% dups_ensembl_remove &
        marker_genes$gene_info$list == 'KM'] <- 'KM+BM'

## Drop the redundant entries
dups_drop <-
    which(
        marker_genes$gene_info$ensembl %in% dups_ensembl_remove &
            marker_genes$gene_info$list == 'BM'
    )
marker_genes$gene_info[dups_drop, ]
#     symbol list     m gene_name         ensembl               models
# 84    RELN   BM  8638      RELN ENSG00000189056               Layer1
# 93    RorB   BM 10825      RORB ENSG00000198963               Layer4
# 101   cux2   BM 14444      CUX2 ENSG00000111249 Layer2+Layer3+Layer4
# 137  Foxp2   BM  8689     FOXP2 ENSG00000128573
marker_genes <- lapply(marker_genes, function(x) {
    x[-dups_drop, ]
})


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


## Process each marker gene (use the full gene matrix to get global stats)
## takes about 3 min to run
layer_idx <- splitit(sce_layer$layer_guess)
eb_markers_list <-
    lapply(seq_len(nrow(marker_genes$gene_info)), function(i) {
        res <- rep(0, ncol(sce_layer))
        layers <-
            colnames(marker_genes$layer_mat)[marker_genes$layer_mat[i, ] > 0]
        res[unlist(layer_idx[layers])] <- 1
        m <- model.matrix( ~ res)
        eBayes(
            lmFit(
                mat,
                design = m,
                block = sce_layer$subject_position,
                correlation = corfit$consensus.correlation
            )
        )
    })


pvals_markers <- sapply(eb_markers_list, function(x) {
    x$p.value[, 2, drop = FALSE]
})

fdr_markers <- apply(pvals_markers, 2, p.adjust, 'fdr')


data.frame(
    'FDRsig' = colSums(apply(pvals0_contrasts, 2, p.adjust, 'fdr') < 0.05),
    'Pval10-6sig' = colSums(pvals0_contrasts < 1e-6),
    'Pval10-8sig' = colSums(pvals0_contrasts < 1e-8)
)

colSums(
    data.frame(
        'FDRsig' = p.adjust(pvals_markers, 'fdr') < 0.05,
        'Pval10-6sig' = pvals_markers < 1e-6,
        'Pval10-8sig' = pvals_markers < 1e-8
    )
)
# FDRsig Pval10.6sig Pval10.8sig
#     87          50          36
length(pvals_markers)
# [1] 130

tstats_markers <- sapply(eb_markers_list, function(x) {
    x$t[, 2, drop = FALSE]
})


## Create a summary
marker_genes_summary <- cbind(marker_genes[[1]], marker_genes[[2]])

## Function to extract the relevant info
marker_extract <- function(m) {
    sapply(seq_along(marker_genes_summary$m), function(i) {
        m[marker_genes_summary$m[i], i]
    })
}

marker_genes_summary$p_value <- marker_extract(pvals_markers)
marker_genes_summary$fdr <- marker_extract(fdr_markers)
marker_genes_summary$t_stat <- marker_extract(tstats_markers)
marker_genes_summary$log_fc <-
    marker_extract(sapply(eb_markers_list, function(x) {
        x$coefficients[, 2, drop = FALSE]
    }))

## Get the ranks in decreasing order, so it's:
# x <- c(-3, 4, 5, 1)
# rank(-x)
## Thus -log10(pval) * -1 = log10(pval)
marker_genes_summary$rank <-
    marker_extract(apply(log10(pvals_markers), 2, rank))

## Or I could have done apply(pvals_markers, 2, rank) and gotten the same result :P
stopifnot(identical(marker_genes_summary$rank,
    marker_extract(apply(
        pvals_markers, 2, rank
    ))))
## And it's the same as the abs t-stats ranks
stopifnot(identical(marker_genes_summary$rank,
    marker_extract(apply(
        -abs(tstats_markers), 2, rank
    ))))


## t-stats should ideally be positive, not all are
table(sign(marker_genes_summary$t_stat))
# -1   1
#  21 105

## The ranks are much better for some results when the t-stat is positive
tapply(marker_genes_summary$rank,
    sign(marker_genes_summary$t_stat),
    summary)
# $`-1`
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    1027    4155    8605    9090   12397   21827
#
# $`1`
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#       1      31     318    3253    2753   21847

## 33 have p-value < 1e-8 for the positive t-stats
tapply(marker_genes_summary$p_value, sign(marker_genes_summary$t_stat), function(x)
    table(x < 1e-8))
# $`-1`
#
# FALSE
#    21
#
# $`1`
#
# FALSE  TRUE
#    72    33

## 47 have p-value < 1e-6 for the positive t-stats
tapply(marker_genes_summary$p_value, sign(marker_genes_summary$t_stat), function(x)
    table(x < 1e-6))
# $`-1`
#
# FALSE
#    21
#
# $`1`
#
# FALSE  TRUE
#    58    47

## 75 have FDR < 0.05 for the positive t-stats
tapply(marker_genes_summary$fdr, sign(marker_genes_summary$t_stat), function(x)
    table(x < 0.05))
# $`-1`
#
# FALSE
#    21
#
# $`1`
#
# FALSE  TRUE
#    30    75


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

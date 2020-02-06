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
    species = c(genes_km_raw$Species, rep('MH', nrow(genes_bm_raw))),
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
#     symbol list species     m gene_name         ensembl               models
# 84    RELN   BM      MH  8638      RELN ENSG00000189056               Layer1
# 93    RorB   BM      MH 10825      RORB ENSG00000198963               Layer4
# 101   cux2   BM      MH 14444      CUX2 ENSG00000111249 Layer2+Layer3+Layer4
# 137  Foxp2   BM      MH  8689     FOXP2 ENSG00000128573               Layer6
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

## What is the estimated block correlation?
corfit$consensus.correlation
# [1] 0.07771005

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

## Extract p-values, FDR, t-stat, logFC
pvals_markers <- sapply(eb_markers_list, function(x) {
    x$p.value[, 2, drop = FALSE]
})
fdr_markers <- apply(pvals_markers, 2, p.adjust, 'fdr')
tstats_markers <- sapply(eb_markers_list, function(x) {
    x$t[, 2, drop = FALSE]
})
logFC_markers <- sapply(eb_markers_list, function(x) {
    x$coefficients[, 2, drop = FALSE]
})


## Create a summary table
marker_genes_summary <- cbind(marker_genes[[1]], marker_genes[[2]])

## Function to extract the relevant info
marker_extract <- function(m, j = marker_genes_summary$m) {
    sapply(seq_along(j), function(i) {
        m[j[i], i]
    })
}

marker_genes_summary$p_value <- marker_extract(pvals_markers)
marker_genes_summary$fdr <- marker_extract(fdr_markers)
marker_genes_summary$t_stat <- marker_extract(tstats_markers)
marker_genes_summary$log_fc <- marker_extract(logFC_markers)

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
addmargins(table(sign(marker_genes_summary$t_stat), marker_genes_summary$species))
#       H  MH Sum
# -1    5  16  21
# 1    17  88 105
# Sum  22 104 126

## The ranks are much better for some results when the t-stat is positive
tapply(marker_genes_summary$rank,
    paste(sign(marker_genes_summary$t_stat), marker_genes_summary$species),
    summary)
# $`-1 H`
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    1792    1977    4155    6862    8980   17405
#
# $`-1 MH`
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    1027    5618    9434    9786   12660   21827
#
# $`1 H`
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#       6      29     738    2039    1883   12092
#
# $`1 MH`
#     Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
#     1.00    46.25   314.50  3487.25  3401.00 21847.00

## 33 have p-value < 1e-8 for the positive t-stats
addmargins(do.call(rbind, tapply(marker_genes_summary$p_value,
    paste(sign(marker_genes_summary$t_stat), marker_genes_summary$species), function(x)
    table(factor(x < 1e-8, levels = c('FALSE', 'TRUE'))))))
#       FALSE TRUE Sum
# -1 H      5    0   5
# -1 MH    16    0  16
# 1 H      13    4  17
# 1 MH     59   29  88
# Sum      93   33 126

## 47 have p-value < 1e-6 for the positive t-stats
addmargins(do.call(rbind, tapply(marker_genes_summary$p_value,
    paste(sign(marker_genes_summary$t_stat), marker_genes_summary$species), function(x)
    table(factor(x < 1e-6, levels = c('FALSE', 'TRUE'))))))
#       FALSE TRUE Sum
# -1 H      5    0   5
# -1 MH    16    0  16
# 1 H       9    8  17
# 1 MH     49   39  88
# Sum      79   47 126

## 75 have FDR < 0.05 for the positive t-stats
addmargins(do.call(rbind, tapply(marker_genes_summary$fdr,
    paste(sign(marker_genes_summary$t_stat), marker_genes_summary$species), function(x)
    table(factor(x < 0.05, levels = c('FALSE', 'TRUE'))))))
#       FALSE TRUE Sum
# -1 H      5    0   5
# -1 MH    16    0  16
# 1 H       5   12  17
# 1 MH     25   63  88
# Sum      51   75 126

## Add mean expr
marker_genes_summary$mean_expr <- rowMeans(mat[marker_genes_summary$m, ])


## Find the best gene for each model
marker_genes_summary$best_gene_m <- apply(log10(pvals_markers), 2, function(x) { which(rank(x) == 1)})
marker_genes_summary$best_gene_name <- rowData(sce_layer)$gene_name[marker_genes_summary$best_gene_m]
marker_genes_summary$best_gene_ensembl <- rownames(sce_layer)[marker_genes_summary$best_gene_m]
marker_genes_summary$best_gene_p_value <- marker_extract(pvals_markers, marker_genes_summary$best_gene_m)
marker_genes_summary$best_gene_fdr <- marker_extract(fdr_markers, marker_genes_summary$best_gene_m)
marker_genes_summary$best_gene_t_stat <- marker_extract(tstats_markers, marker_genes_summary$best_gene_m)
marker_genes_summary$best_gene_log_fc <- marker_extract(logFC_markers, marker_genes_summary$best_gene_m)
marker_genes_summary$best_gene_mean_expr <- rowMeans(mat[marker_genes_summary$best_gene_m, ])

## Re-order by p-value
marker_genes <- lapply(marker_genes, function(x) {
    x[order(marker_genes_summary$p_value), ]
})
marker_genes_summary <- marker_genes_summary[order(marker_genes_summary$p_value), ]

## Save for later
save(eb_markers_list, file = 'rda/eb_markers_list.Rdata')
save(marker_genes, marker_genes_summary, file = 'rda/marker_genes_summary.Rdata')

write.csv(marker_genes_summary,
    file = 'marker_genes_summary.csv',
    quote = FALSE,
    row.names = FALSE)


## Make boxplots
layer_guess_reordered <-  factor(sce_layer$layer_guess, levels = c(paste0('Layer', 1:6), 'WM'))
pdf('pdf/markers_gene_optimal_boxplots.pdf', useDingbats = FALSE)
set.seed(20200206)
for (i in seq_len(nrow(marker_genes_summary))) {
    # i <- 1
    message(paste(Sys.time(), 'making the plot for', i, 'gene', marker_genes_summary$gene_name[i]))
    curr_layers <- strsplit(marker_genes_summary$models[i], '\\+')[[1]]
    boxplot(
        mat[marker_genes_summary$m[i],] ~ layer_guess_reordered,
        xlab = 'Layer',
        ylab = 'logcounts',
        main = paste(
            marker_genes_summary$gene_name[i],
            marker_genes_summary$ensembl[i],
            paste(ifelse(marker_genes_summary$species[i] == 'H', 'human', 'mouse'), 'marker'),
            paste0('(', marker_genes_summary$list[i], ')'),
            '\n',
            paste(marker_genes_summary$models[i], 'vs others'),
            '\n',
            'tstat',
            formatC(marker_genes_summary$t_stat[i], format = "e", digits = 2),
            'p',
            formatC(marker_genes_summary$p_val[i], format = "e", digits = 2),
            'rank',
            marker_genes_summary$rank[i]
        ),
        outline = FALSE,
        cex = 1.5,
        col = ifelse(levels(layer_guess_reordered) %in% curr_layers, 'skyblue', 'violet')
    )
    points(
        mat[marker_genes_summary$m[i],] ~ jitter(as.integer(layer_guess_reordered)),
        pch = 21,
        bg = ifelse(layer_guess_reordered %in% curr_layers, 'dodgerblue4', 'darkviolet'),
        cex = 1.5
    )
    # legend(
    #     "top",
    #     levels(sce_layer$layer_guess),
    #     col =  Polychrome::palette36.colors(7),
    #     lwd = 2
    # )
}
dev.off()



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

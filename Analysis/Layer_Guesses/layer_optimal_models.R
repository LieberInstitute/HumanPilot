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
library('ggplot2')

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
        rowSums(genes_query_mat) > 0,],
    layer_mat = genes_query_mat[!is.na(genes_query$m) &
            rowSums(genes_query_mat) > 0,])

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
marker_genes$gene_info[dups_drop,]
#     symbol list species     m gene_name         ensembl               models
# 84    RELN   BM      MH  8638      RELN ENSG00000189056               Layer1
# 93    RorB   BM      MH 10825      RORB ENSG00000198963               Layer4
# 101   cux2   BM      MH 14444      CUX2 ENSG00000111249 Layer2+Layer3+Layer4
# 137  Foxp2   BM      MH  8689     FOXP2 ENSG00000128573               Layer6
marker_genes <- lapply(marker_genes, function(x) {
    x[-dups_drop,]
})


##### Compute the duplicate correlation globally
## Extract the data
mat <- assays(sce_layer)$logcounts

## Build a group model
mod <- with(colData(sce_layer), model.matrix(~ 0 + layer_guess))
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
            colnames(marker_genes$layer_mat)[marker_genes$layer_mat[i,] > 0]
        res[unlist(layer_idx[layers])] <- 1
        m <- model.matrix(~ res)
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
addmargins(table(
    sign(marker_genes_summary$t_stat),
    marker_genes_summary$species
))
#       H  MH Sum
# -1    5  16  21
# 1    17  88 105
# Sum  22 104 126

## The ranks are much better for some results when the t-stat is positive
tapply(marker_genes_summary$rank,
    paste(
        sign(marker_genes_summary$t_stat),
        marker_genes_summary$species
    ),
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
addmargins(do.call(
    rbind,
    tapply(marker_genes_summary$p_value,
        paste(
            sign(marker_genes_summary$t_stat),
            marker_genes_summary$species
        ), function(x)
            table(factor(x < 1e-8, levels = c('FALSE', 'TRUE'))))
))
#       FALSE TRUE Sum
# -1 H      5    0   5
# -1 MH    16    0  16
# 1 H      13    4  17
# 1 MH     59   29  88
# Sum      93   33 126

## 47 have p-value < 1e-6 for the positive t-stats
addmargins(do.call(
    rbind,
    tapply(marker_genes_summary$p_value,
        paste(
            sign(marker_genes_summary$t_stat),
            marker_genes_summary$species
        ), function(x)
            table(factor(x < 1e-6, levels = c('FALSE', 'TRUE'))))
))
#       FALSE TRUE Sum
# -1 H      5    0   5
# -1 MH    16    0  16
# 1 H       9    8  17
# 1 MH     49   39  88
# Sum      79   47 126

## 75 have FDR < 0.05 for the positive t-stats
addmargins(do.call(
    rbind,
    tapply(marker_genes_summary$fdr,
        paste(
            sign(marker_genes_summary$t_stat),
            marker_genes_summary$species
        ), function(x)
            table(factor(x < 0.05, levels = c('FALSE', 'TRUE'))))
))
#       FALSE TRUE Sum
# -1 H      5    0   5
# -1 MH    16    0  16
# 1 H       5   12  17
# 1 MH     25   63  88
# Sum      51   75 126

## Add mean expr
marker_genes_summary$mean_expr <-
    rowMeans(mat[marker_genes_summary$m,])


## Find the best gene for each model
marker_genes_summary$best_gene_m <-
    apply(log10(pvals_markers), 2, function(x) {
        which(rank(x) == 1)
    })
marker_genes_summary$best_gene_name <-
    rowData(sce_layer)$gene_name[marker_genes_summary$best_gene_m]
marker_genes_summary$best_gene_ensembl <-
    rownames(sce_layer)[marker_genes_summary$best_gene_m]
marker_genes_summary$best_gene_p_value <-
    marker_extract(pvals_markers, marker_genes_summary$best_gene_m)
marker_genes_summary$best_gene_fdr <-
    marker_extract(fdr_markers, marker_genes_summary$best_gene_m)
marker_genes_summary$best_gene_t_stat <-
    marker_extract(tstats_markers, marker_genes_summary$best_gene_m)
marker_genes_summary$best_gene_log_fc <-
    marker_extract(logFC_markers, marker_genes_summary$best_gene_m)
marker_genes_summary$best_gene_mean_expr <-
    rowMeans(mat[marker_genes_summary$best_gene_m,])

## Re-order by p-value
marker_genes <- lapply(marker_genes, function(x) {
    x[order(marker_genes_summary$p_value),]
})
marker_genes_summary <-
    marker_genes_summary[order(marker_genes_summary$p_value),]

## Save for later
save(eb_markers_list, file = 'rda/eb_markers_list.Rdata')
save(marker_genes, marker_genes_summary, file = 'rda/marker_genes_summary.Rdata')

write.csv(
    marker_genes_summary,
    file = 'marker_genes_summary.csv',
    quote = FALSE,
    row.names = FALSE
)


## Make boxplots
layer_guess_reordered <-
    factor(sce_layer$layer_guess, levels = c(paste0('Layer', 1:6), 'WM'))
pdf('pdf/markers_gene_optimal_boxplots.pdf', useDingbats = FALSE)
set.seed(20200206)
for (i in seq_len(nrow(marker_genes_summary))) {
    # i <- 1
    message(
        paste(
            Sys.time(),
            'making the plot for',
            i,
            'gene',
            marker_genes_summary$gene_name[i]
        )
    )
    curr_layers <-
        strsplit(marker_genes_summary$models[i], '\\+')[[1]]
    boxplot(
        mat[marker_genes_summary$m[i], ] ~ layer_guess_reordered,
        xlab = 'Layer',
        ylab = 'logcounts',
        main = paste(
            marker_genes_summary$gene_name[i],
            marker_genes_summary$ensembl[i],
            paste(
                ifelse(marker_genes_summary$species[i] == 'H', 'human', 'mouse'),
                'marker'
            ),
            paste0('(', marker_genes_summary$list[i], ')'),
            '\n',
            paste(marker_genes_summary$models[i], 'vs others'),
            '\n',
            'tstat',
            formatC(
                marker_genes_summary$t_stat[i],
                format = "e",
                digits = 2
            ),
            'p',
            formatC(
                marker_genes_summary$p_val[i],
                format = "e",
                digits = 2
            ),
            'rank',
            marker_genes_summary$rank[i]
        ),
        outline = FALSE,
        cex = 1.5,
        col = ifelse(
            levels(layer_guess_reordered) %in% curr_layers,
            'skyblue',
            'violet'
        )
    )
    points(
        mat[marker_genes_summary$m[i], ] ~ jitter(as.integer(layer_guess_reordered)),
        pch = 21,
        bg = ifelse(
            layer_guess_reordered %in% curr_layers,
            'dodgerblue4',
            'darkviolet'
        ),
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


marker_genes_summary$species <-
    ifelse(marker_genes_summary$species == 'MH', 'mouse', 'human')

pdf('pdf/markers_gene_optimal_scatter.pdf', useDingbats = FALSE)
ggplot(marker_genes_summary,
    aes(
        x = -log10(best_gene_p_value),
        y = -log10(p_value),
        shape = species,
        color = rank < 100
    )) + geom_point() + theme_bw(base_size = 20) + ylab('-log10 p-value') + xlab('-log10 p-value for the best gene')
ggplot(marker_genes_summary,
    aes(
        x = -log10(best_gene_fdr),
        y = -log10(fdr),
        shape = species,
        color = rank < 100
    )) + geom_point() + theme_bw(base_size = 20) + ylab('-log10 FDR') + xlab('-log10 FDR for the best gene')
ggplot(
    marker_genes_summary,
    aes(
        x = best_gene_t_stat,
        y = t_stat,
        shape = species,
        color = rank < 100
    )
) + geom_point() + coord_fixed() + theme_bw(base_size = 20) + ylab('t-statistic') + xlab('t-statistic for the best gene')
ggplot(
    marker_genes_summary,
    aes(
        x = best_gene_log_fc,
        y = log_fc,
        shape = species,
        color = rank < 100
    )
) + geom_point() + coord_fixed() + theme_bw(base_size = 20) + ylab('log fold change') + xlab('log fold change for the best gene')
dev.off()


with(marker_genes_summary, cor(-log10(best_gene_p_value), -log10(p_value)))
# [1] 0.4224215
with(marker_genes_summary, cor(-log10(best_gene_fdr), -log10(fdr)))
# [1] 0.4719611
with(marker_genes_summary, cor(best_gene_t_stat, t_stat))
# [1] 0.144047
with(marker_genes_summary, cor(best_gene_log_fc, log_fc))
# [1] 0.1767757

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 3.6.1 Patched (2019-10-31 r77350)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2020-02-06
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  assertthat             0.2.1     2019-03-21 [2] CRAN (R 3.6.1)
#  backports              1.1.5     2019-10-02 [1] CRAN (R 3.6.1)
#  beeswarm               0.2.3     2016-04-25 [2] CRAN (R 3.6.1)
#  Biobase              * 2.46.0    2019-10-29 [2] Bioconductor
#  BiocGenerics         * 0.32.0    2019-10-29 [1] Bioconductor
#  BiocNeighbors          1.4.1     2019-11-01 [2] Bioconductor
#  BiocParallel         * 1.20.1    2019-12-21 [1] Bioconductor
#  BiocSingular           1.2.0     2019-10-29 [2] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.6.1)
#  cellranger             1.1.0     2016-07-27 [1] CRAN (R 3.6.1)
#  cli                    2.0.0     2019-12-09 [1] CRAN (R 3.6.1)
#  cluster              * 2.1.0     2019-06-19 [3] CRAN (R 3.6.1)
#  colorout             * 1.2-2     2019-10-31 [1] Github (jalvesaq/colorout@641ed38)
#  colorspace             1.4-1     2019-03-18 [2] CRAN (R 3.6.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.6.1)
#  DelayedArray         * 0.12.0    2019-10-29 [2] Bioconductor
#  DelayedMatrixStats     1.8.0     2019-10-29 [2] Bioconductor
#  digest                 0.6.23    2019-11-23 [1] CRAN (R 3.6.1)
#  dplyr                  0.8.3     2019-07-04 [1] CRAN (R 3.6.1)
#  dqrng                  0.2.1     2019-05-17 [1] CRAN (R 3.6.1)
#  edgeR                  3.28.0    2019-10-29 [1] Bioconductor
#  fansi                  0.4.0     2018-10-05 [1] CRAN (R 3.6.1)
#  GenomeInfoDb         * 1.22.0    2019-10-29 [1] Bioconductor
#  GenomeInfoDbData       1.2.2     2019-10-28 [2] Bioconductor
#  GenomicRanges        * 1.38.0    2019-10-29 [1] Bioconductor
#  ggbeeswarm             0.6.0     2017-08-07 [2] CRAN (R 3.6.1)
#  ggplot2              * 3.2.1     2019-08-10 [1] CRAN (R 3.6.1)
#  glue                   1.3.1     2019-03-12 [1] CRAN (R 3.6.1)
#  googledrive            1.0.0     2019-08-19 [1] CRAN (R 3.6.1)
#  gridExtra              2.3       2017-09-09 [2] CRAN (R 3.6.1)
#  gtable                 0.3.0     2019-03-25 [2] CRAN (R 3.6.1)
#  here                 * 0.1       2017-05-28 [1] CRAN (R 3.6.1)
#  htmltools              0.4.0     2019-10-04 [1] CRAN (R 3.6.1)
#  htmlwidgets            1.5.1     2019-10-08 [1] CRAN (R 3.6.1)
#  httpuv                 1.5.2     2019-09-11 [1] CRAN (R 3.6.1)
#  igraph                 1.2.4.1   2019-04-22 [2] CRAN (R 3.6.1)
#  IRanges              * 2.20.1    2019-11-20 [1] Bioconductor
#  irlba                  2.3.3     2019-02-05 [1] CRAN (R 3.6.1)
#  jaffelab             * 0.99.29   2019-11-04 [1] Github (LieberInstitute/jaffelab@a7d87cb)
#  jsonlite               1.6       2018-12-07 [2] CRAN (R 3.6.1)
#  labeling               0.3       2014-08-23 [2] CRAN (R 3.6.1)
#  later                  1.0.0     2019-10-04 [1] CRAN (R 3.6.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.6.1)
#  lazyeval               0.2.2     2019-03-15 [2] CRAN (R 3.6.1)
#  limma                * 3.42.0    2019-10-29 [1] Bioconductor
#  locfit                 1.5-9.1   2013-04-20 [2] CRAN (R 3.6.1)
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.6.1)
#  Matrix                 1.2-17    2019-03-22 [3] CRAN (R 3.6.1)
#  matrixStats          * 0.55.0    2019-09-07 [1] CRAN (R 3.6.1)
#  mime                   0.8       2019-12-19 [1] CRAN (R 3.6.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.6.1)
#  pheatmap             * 1.0.12    2019-01-04 [2] CRAN (R 3.6.1)
#  pillar                 1.4.3     2019-12-20 [1] CRAN (R 3.6.1)
#  pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 3.6.1)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.6.1)
#  Polychrome           * 1.2.3     2019-08-01 [1] CRAN (R 3.6.1)
#  promises               1.1.0     2019-10-04 [1] CRAN (R 3.6.1)
#  purrr                  0.3.3     2019-10-18 [2] CRAN (R 3.6.1)
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.6.1)
#  rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 3.6.1)
#  RColorBrewer           1.1-2     2014-12-07 [2] CRAN (R 3.6.1)
#  Rcpp                   1.0.3     2019-11-08 [1] CRAN (R 3.6.1)
#  RCurl                  1.95-4.12 2019-03-04 [2] CRAN (R 3.6.1)
#  readxl               * 1.3.1     2019-03-13 [2] CRAN (R 3.6.1)
#  rlang                  0.4.2     2019-11-23 [1] CRAN (R 3.6.1)
#  rmote                * 0.3.4     2019-10-31 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot              1.3-2     2018-01-03 [2] CRAN (R 3.6.1)
#  rsvd                   1.0.2     2019-07-29 [1] CRAN (R 3.6.1)
#  S4Vectors            * 0.24.1    2019-12-01 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.6.1)
#  scater               * 1.14.3    2019-11-07 [2] Bioconductor
#  scatterplot3d          0.3-41    2018-03-14 [1] CRAN (R 3.6.1)
#  scran                * 1.14.5    2019-11-19 [1] Bioconductor
#  segmented              1.0-0     2019-06-17 [2] CRAN (R 3.6.1)
#  servr                  0.15      2019-08-07 [1] CRAN (R 3.6.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.6.1)
#  SingleCellExperiment * 1.8.0     2019-10-29 [2] Bioconductor
#  statmod                1.4.32    2019-05-29 [2] CRAN (R 3.6.1)
#  SummarizedExperiment * 1.16.1    2019-12-19 [1] Bioconductor
#  tibble                 2.1.3     2019-06-06 [1] CRAN (R 3.6.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.6.1)
#  utf8                   1.1.4     2018-05-24 [1] CRAN (R 3.6.1)
#  vctrs                  0.2.1     2019-12-17 [1] CRAN (R 3.6.1)
#  vipor                  0.4.5     2017-03-22 [2] CRAN (R 3.6.1)
#  viridis                0.5.1     2018-03-29 [2] CRAN (R 3.6.1)
#  viridisLite            0.3.0     2018-02-01 [2] CRAN (R 3.6.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.6.1)
#  xfun                   0.11      2019-11-12 [1] CRAN (R 3.6.1)
#  XVector                0.26.0    2019-10-29 [1] Bioconductor
#  zeallot                0.1.0     2018-01-28 [1] CRAN (R 3.6.1)
#  zlibbioc               1.32.0    2019-10-29 [2] Bioconductor
#
# [1] /users/lcollado/R/3.6.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/library

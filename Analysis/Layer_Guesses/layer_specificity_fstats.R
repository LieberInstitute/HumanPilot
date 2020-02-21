library('SingleCellExperiment')
library('here')
library('readxl')
library('limma')
library('sessioninfo')

dir.create('pdf', showWarnings = FALSE)
dir.create('rda', showWarnings = FALSE)

## Load data
load(here(
    'Analysis',
    'Human_DLPFC_Visium_processedData_sce_scran.Rdata'
))

## Functions derived from this script, to make it easier to resume the work
sce_layer_file <-
    here('Analysis', 'Layer_Guesses', 'rda', 'sce_layer.Rdata')
if (file.exists(sce_layer_file))
    load(sce_layer_file, verbose = TRUE)
source(here('Analysis', 'Layer_Guesses', 'layer_specificity_functions.R'))


## From layer_specificity.R
fit_f_model <- function(sce) {
    message(paste(Sys.time(), 'starting the model run'))

    ## Extract the data
    mat <- assays(sce)$logcounts

    ## For dropping un-used levels
    sce$layer_guess <- factor(sce$layer_guess)

    ## Build a group model
    mod <- with(colData(sce), model.matrix( ~ layer_guess))
    colnames(mod) <- gsub('layer_guess', '', colnames(mod))
    ## Takes like 2 min to run
    corfit <-
        duplicateCorrelation(mat, mod, block = sce$subject_position)
    message(paste(Sys.time(), 'correlation:', corfit$consensus.correlation))
    fit <-
        lmFit(
            mat,
            design = mod,
            block = sce$subject_position,
            correlation = corfit$consensus.correlation
        )
    eb <- eBayes(fit)
    return(eb)
}

ebF_list <-
    lapply(list('full' = sce_layer, 'noWM' = sce_layer[, sce_layer$layer_guess != 'WM']), fit_f_model)
# 2020-02-21 10:13:31 starting the model run
# 2020-02-21 10:14:27 correlation: 0.0777100513762501
# 2020-02-21 10:14:27 starting the model run
# 2020-02-21 10:15:19 correlation: 0.0959024426789848


## Extract F-statistics
f_stats <- do.call(cbind, lapply(names(ebF_list), function(i) {
    x <- ebF_list[[i]]
    top <-
        topTable(
            x,
            coef = 2:ncol(x$coefficients),
            sort.by = 'none',
            number = length(x$F)
        )
    # identical(p.adjust(top$P.Value, 'fdr'), top$adj.P.Val)
    res <- data.frame(
        'f' = top$F,
        'p_value' = top$P.Value,
        'fdr' = top$adj.P.Val,
        'AveExpr' = top$AveExpr,
        stringsAsFactors = FALSE
    )
    colnames(res) <- paste0(i, '_', colnames(res))
    return(res)
}))
f_stats$ensembl <- rownames(sce_layer)
f_stats$gene <- rowData(sce_layer)$gene_name
rownames(f_stats) <- NULL

head(f_stats)
#      full_f full_p_value   full_fdr full_AveExpr    noWM_f noWM_p_value
# 1 0.7055034   0.64613823 0.71662018    0.1422407 0.6753082   0.64375784
# 2 2.5500608   0.02705036 0.05586955    0.8383094 1.6289143   0.16611181
# 3 2.2086804   0.05191480 0.09724935    3.5454133 2.1242007   0.07466448
# 4 0.9785169   0.44622269 0.54759570    0.1048252 0.7183197   0.61218520
# 5 1.6674804   0.14158125 0.22362788    1.2135801 1.1911356   0.32429532
# 6 1.3364877   0.25244796 0.35390894    1.3447335 0.5876680   0.70932745
#    noWM_fdr noWM_AveExpr         ensembl        gene
# 1 0.7283289    0.1565079 ENSG00000243485 MIR1302-2HG
# 2 0.2774245    0.9362990 ENSG00000238009  AL627309.1
# 3 0.1498954    3.6730061 ENSG00000237491  AL669831.5
# 4 0.7040224    0.1244799 ENSG00000177757      FAM87B
# 5 0.4526797    1.1121321 ENSG00000225880   LINC00115
# 6 0.7800006    1.4543898 ENSG00000230368      FAM41C

pdf('pdf/layer_specificity_full_vs_noWM.pdf', useDingbats = FALSE)
with(
    f_stats,
    plot(
        full_f,
        noWM_f,
        xlab = 'F-stats (WM + L1 through L6)',
        ylab = 'F-stats (L1 through L6 only)',
        pch = 21,
        bg = add.alpha('black', 1 / 5),
        col = add.alpha('black', 1 / 5),
        cex = 0.8,
        cex.lab = 1.3
    )
)
abline(a = 0, b = 1, col = 'red')

with(
    f_stats,
    plot(
        -log10(full_p_value),
        -log10(noWM_p_value),
        xlab = '-log10 p-value (WM + L1 through L6)',
        ylab = '-log10 p-value (L1 through L6 only)',
        pch = 21,
        bg = add.alpha('black', 1 / 5),
        col = add.alpha('black', 1 / 5),
        cex = 0.8,
        cex.lab = 1.3
    )
)
abline(a = 0, b = 1, col = 'red')
dev.off()

## Move the ones that we made originally with the 0+layer_guess model
# system('mv rda/ebF_list.Rdata rda/ebF_list_model_0_plus_layer_guess.Rdata')
## Save for later
save(f_stats, ebF_list, file = 'rda/ebF_list.Rdata')

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
#  date     2020-02-21
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version  date       lib source
#  assertthat             0.2.1    2019-03-21 [2] CRAN (R 3.6.1)
#  backports              1.1.5    2019-10-02 [1] CRAN (R 3.6.1)
#  Biobase              * 2.46.0   2019-10-29 [2] Bioconductor
#  BiocGenerics         * 0.32.0   2019-10-29 [1] Bioconductor
#  BiocParallel         * 1.20.1   2019-12-21 [1] Bioconductor
#  bitops                 1.0-6    2013-08-17 [2] CRAN (R 3.6.1)
#  cellranger             1.1.0    2016-07-27 [1] CRAN (R 3.6.1)
#  cli                    2.0.0    2019-12-09 [1] CRAN (R 3.6.1)
#  colorout             * 1.2-2    2019-10-31 [1] Github (jalvesaq/colorout@641ed38)
#  colorspace             1.4-1    2019-03-18 [2] CRAN (R 3.6.1)
#  crayon                 1.3.4    2017-09-16 [1] CRAN (R 3.6.1)
#  DelayedArray         * 0.12.2   2020-01-06 [2] Bioconductor
#  digest                 0.6.23   2019-11-23 [1] CRAN (R 3.6.1)
#  fansi                  0.4.0    2018-10-05 [1] CRAN (R 3.6.1)
#  GenomeInfoDb         * 1.22.0   2019-10-29 [1] Bioconductor
#  GenomeInfoDbData       1.2.2    2019-10-28 [2] Bioconductor
#  GenomicRanges        * 1.38.0   2019-10-29 [1] Bioconductor
#  ggplot2                3.2.1    2019-08-10 [1] CRAN (R 3.6.1)
#  glue                   1.3.1    2019-03-12 [1] CRAN (R 3.6.1)
#  gtable                 0.3.0    2019-03-25 [2] CRAN (R 3.6.1)
#  here                 * 0.1      2017-05-28 [1] CRAN (R 3.6.1)
#  htmltools              0.4.0    2019-10-04 [1] CRAN (R 3.6.1)
#  htmlwidgets            1.5.1    2019-10-08 [1] CRAN (R 3.6.1)
#  httpuv                 1.5.2    2019-09-11 [1] CRAN (R 3.6.1)
#  IRanges              * 2.20.1   2019-11-20 [1] Bioconductor
#  jsonlite               1.6.1    2020-02-02 [2] CRAN (R 3.6.1)
#  later                  1.0.0    2019-10-04 [1] CRAN (R 3.6.1)
#  lattice                0.20-38  2018-11-04 [3] CRAN (R 3.6.1)
#  lazyeval               0.2.2    2019-03-15 [2] CRAN (R 3.6.1)
#  lifecycle              0.1.0    2019-08-01 [1] CRAN (R 3.6.1)
#  limma                * 3.42.0   2019-10-29 [1] Bioconductor
#  magrittr               1.5      2014-11-22 [1] CRAN (R 3.6.1)
#  Matrix                 1.2-17   2019-03-22 [3] CRAN (R 3.6.1)
#  matrixStats          * 0.55.0   2019-09-07 [1] CRAN (R 3.6.1)
#  munsell                0.5.0    2018-06-12 [2] CRAN (R 3.6.1)
#  pillar                 1.4.3    2019-12-20 [1] CRAN (R 3.6.1)
#  pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 3.6.1)
#  png                    0.1-7    2013-12-03 [2] CRAN (R 3.6.1)
#  promises               1.1.0    2019-10-04 [1] CRAN (R 3.6.1)
#  R6                     2.4.1    2019-11-12 [2] CRAN (R 3.6.1)
#  RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 3.6.1)
#  Rcpp                   1.0.3    2019-11-08 [1] CRAN (R 3.6.1)
#  RCurl                  1.98-1.1 2020-01-19 [2] CRAN (R 3.6.1)
#  readxl               * 1.3.1    2019-03-13 [2] CRAN (R 3.6.1)
#  rlang                  0.4.2    2019-11-23 [1] CRAN (R 3.6.1)
#  rmote                * 0.3.4    2019-10-31 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot              1.3-2    2018-01-03 [2] CRAN (R 3.6.1)
#  S4Vectors            * 0.24.1   2019-12-01 [1] Bioconductor
#  scales                 1.1.0    2019-11-18 [2] CRAN (R 3.6.1)
#  servr                  0.15     2019-08-07 [1] CRAN (R 3.6.1)
#  sessioninfo          * 1.1.1    2018-11-05 [1] CRAN (R 3.6.1)
#  SingleCellExperiment * 1.8.0    2019-10-29 [2] Bioconductor
#  statmod                1.4.34   2020-02-17 [2] CRAN (R 3.6.1)
#  SummarizedExperiment * 1.16.1   2019-12-19 [1] Bioconductor
#  tibble                 2.1.3    2019-06-06 [1] CRAN (R 3.6.1)
#  withr                  2.1.2    2018-03-15 [2] CRAN (R 3.6.1)
#  xfun                   0.11     2019-11-12 [1] CRAN (R 3.6.1)
#  XVector                0.26.0   2019-10-29 [1] Bioconductor
#  zlibbioc               1.32.0   2019-10-29 [2] Bioconductor
#
# [1] /users/lcollado/R/3.6.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/library

library('SingleCellExperiment')
library('here')
library('readxl')
library('limma')
library('grid')
library('gridExtra')
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

## For plotting
source(here('Analysis', 'spatialLIBD_global_plot_code.R'))
genes <- paste0(rowData(sce)$gene_name, '; ', rowData(sce)$gene_id)


## mean Xk unique molecular indices (UMIs) and mean Xk genes per spot
summary(sce$sum_umi)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   17    2035    3034    3462    4407   20600

summary(sce$sum_gene)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   16    1178    1631    1734    2176    6035


## chrM genes
ix_mito <- grep("^MT-", rowData(sce)$gene_name)
rowData(sce)$gene_name[ix_mito]
# [1] "MT-ND1"  "MT-ND2"  "MT-CO1"  "MT-CO2"  "MT-ATP8" "MT-ATP6" "MT-CO3"
# [8] "MT-ND3"  "MT-ND4L" "MT-ND4"  "MT-ND5"  "MT-ND6"  "MT-CYB"

## Should save this on the sce object later
expr_total <- colSums(assays(sce)$counts)
## Actually, we already had this
identical(sce$sum_umi, expr_total)
# [1] TRUE
expr_chrM <- colSums(assays(sce)$counts[ix_mito,])
expr_chrM_ratio <- expr_chrM / expr_total
summary(expr_chrM_ratio)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.04853 0.15465 0.18442 0.18554 0.21521 0.44156


## Visualize this data at the spot level
## In the future we could customize the colors if we want to
sce$expr_total <- expr_total
sce$expr_chrM <- expr_chrM
sce$expr_chrM_ratio <- expr_chrM_ratio
sce_image_grid_gene(
    sce,
    geneid = 'expr_total',
    spatial = TRUE,
    minCount = 0,
    pdf_file = 'pdf/spot_expr_total.pdf'
)
sce_image_grid_gene(
    sce,
    geneid = 'expr_chrM',
    spatial = TRUE,
    minCount = 0,
    pdf_file = 'pdf/spot_expr_chrM.pdf'
)
sce_image_grid_gene(
    sce,
    geneid = 'expr_chrM_ratio',
    spatial = TRUE,
    minCount = 0,
    pdf_file = 'pdf/spot_expr_chrM_ratio.pdf'
)


## Repeat at the layer-level
# ix_mito_layer <- grep("^MT-", rowData(sce_layer)$gene_name)
# expr_total_layer <- colSums(assays(sce_layer)$counts)
# expr_chrM_layer <-
#     colSums(assays(sce_layer)$counts[ix_mito_layer,])
# expr_chrM_ratio_layer <- expr_chrM_layer / expr_total_layer
# summary(expr_chrM_ratio_layer)
## Err, it's all 0 because we already dropped chrM by this point :P


## chrM_ratio vs number of cells
pdf('pdf/spot_expr_chrM_ratio_by_cell_boxplot.pdf',
    useDingbats = FALSE)

## From https://community.rstudio.com/t/add-regression-line-in-geom-boxplot/9096/2?u=lcolladotor
coefs <- coef(lm(sce$expr_chrM_ratio ~ sce$cell_count))
ggplot(as.data.frame(colData(sce)),
    aes(x = cell_count, group = cell_count, y = expr_chrM_ratio)) +
    geom_boxplot() + ylab('chrM expression ratio') +
    xlab('Number of cells per spot') +
    theme_bw(base_size = 20) + xlim(c(-1, 30)) +
    geom_abline(
        intercept = coefs[1],
        slope = coefs[2],
        colour = 'red',
        linetype = 2
    )
dev.off()


f_chrM_cell <- function(sample = 1) {
    sce <- sce[, sce$sample_name == unique(sce$sample_name)[sample]]
    coefs <- coef(lm(sce$expr_chrM_ratio ~ sce$cell_count))
    ggplot(
        as.data.frame(colData(sce)),
        aes(x = cell_count, group = cell_count, y = expr_chrM_ratio)
    ) +
        geom_boxplot() + ylab('') +
        xlab('') +
        theme_bw(base_size = 20) + xlim(c(-1, 30)) +
        facet_grid( ~ sample_name) +
        geom_abline(
            intercept = coefs[1],
            slope = coefs[2],
            colour = 'red',
            linetype = 2
        )
}

p_list <- cowplot::plot_grid(
    plotlist = lapply(1:12, f_chrM_cell),
    nrow = 3,
    ncol = 4
)
pdf(
    'pdf/spot_expr_chrM_ratio_by_cell_boxplot_per_sample.pdf',
    useDingbats = FALSE,
    width = 4 * 5,
    height = 3 * 5
)
## From https://stackoverflow.com/questions/33114380/centered-x-axis-label-for-muliplot-using-cowplot-package
#create common x and y labels
y.grob <- textGrob(
    "chrM expression ratio",
    gp = gpar(
        fontface = "bold",
        col = "black",
        fontsize = 40
    ),
    rot = 90
)
x.grob <- textGrob("Number of cells per spot",
    gp = gpar(
        fontface = "bold",
        col = "black",
        fontsize = 40
    ))

#add to plot
grid.arrange(arrangeGrob(p_list, left = y.grob, bottom = x.grob))
dev.off()

## Visualize some genes
genes[grep('SNAP25', genes)]
# [1] "SNAP25-AS1; ENSG00000227906"

sce_image_grid_gene(
    sce,
    geneid = "SNAP25; ENSG00000132639",
    spatial = TRUE,
    minCount = 0,
    pdf_file = 'pdf/spot_expr_SNAP25.pdf'
)

genes[grep('MOBP', genes)]
# [1] "MOBP; ENSG00000168314"

sce_image_grid_gene(
    sce,
    geneid = "MOBP; ENSG00000168314",
    spatial = TRUE,
    minCount = 0,
    pdf_file = 'pdf/spot_expr_MOBP.pdf'
)



### Check some stat outputs

## load modeling outputs
load("rda/eb_contrasts.Rdata", verbose = TRUE)
load("rda/eb0_list.Rdata", verbose = TRUE)
load('rda/ebF_list.Rdata', verbose = TRUE)

### Compute these numbers
## with XXXX DE genes (DEGs) across the seven layers (at FDR < 0.05) and
## XXXX genes across the six layers (excluding white matter, at FDR < 0.05)

f_sig <- function(type, cut = 0.05) {
    cbind('n' = addmargins(table(f_stats[[type]] < cut)),
        'ratio' = addmargins(table(f_stats[[type]] < cut)) / nrow(f_stats))
}
f_sig('full_fdr')
#           n     ratio
# FALSE  2433 0.1089517
# TRUE  19898 0.8910483
# Sum   22331 1.0000000

f_sig('noWM_fdr')
#           n     ratio
# FALSE  2563 0.1147732
# TRUE  19768 0.8852268
# Sum   22331 1.0000000

f_sig('full_fdr', 0.001)
#           n    ratio
# FALSE  4220 0.188975
# TRUE  18111 0.811025
# Sum   22331 1.000000

f_sig('noWM_fdr', 0.001)
#           n    ratio
# FALSE  4393 0.196722
# TRUE  17938 0.803278
# Sum   22331 1.000000

## Extract the p-values
pvals0_contrasts <- sapply(eb0_list, function(x) {
    x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts) = rownames(eb_contrasts)
fdrs0_contrasts = apply(pvals0_contrasts, 2, p.adjust, "fdr")

## Extract the t-stats
t0_contrasts <- sapply(eb0_list, function(x) {
    x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts) = rownames(eb_contrasts)

summary(fdrs0_contrasts < 0.05)
#     WM            Layer1          Layer2          Layer3
# Mode :logical   Mode :logical   Mode :logical   Mode :logical
# FALSE:13207     FALSE:19298     FALSE:20769     FALSE:22148
# TRUE :9124      TRUE :3033      TRUE :1562      TRUE :183
#   Layer4          Layer5          Layer6
# Mode :logical   Mode :logical   Mode :logical
# FALSE:21591     FALSE:21688     FALSE:21952
# TRUE :740       TRUE :643       TRUE :379
sort(colSums(fdrs0_contrasts < 0.05))
# Layer3 Layer6 Layer5 Layer4 Layer2 Layer1     WM
#    183    379    643    740   1562   3033   9124

pvals_contrasts <- eb_contrasts$p.value
fdrs_contrasts <- apply(pvals_contrasts, 2, p.adjust, "fdr")
dim(pvals_contrasts)
# [1] 22331    21

summary(fdrs_contrasts < 0.05)
# WM-Layer1       WM-Layer2       WM-Layer3       WM-Layer4
#  Mode :logical   Mode :logical   Mode :logical   Mode :logical
#  FALSE:16664     FALSE:14339     FALSE:13831     FALSE:13873
#  TRUE :5667      TRUE :7992      TRUE :8500      TRUE :8458
#  WM-Layer5       WM-Layer6       Layer1-Layer2   Layer1-Layer3
#  Mode :logical   Mode :logical   Mode :logical   Mode :logical
#  FALSE:14352     FALSE:15645     FALSE:18645     FALSE:18765
#  TRUE :7979      TRUE :6686      TRUE :3686      TRUE :3566
#  Layer1-Layer4   Layer1-Layer5   Layer1-Layer6   Layer2-Layer3
#  Mode :logical   Mode :logical   Mode :logical   Mode :logical
#  FALSE:17654     FALSE:17693     FALSE:18076     FALSE:21954
#  TRUE :4677      TRUE :4638      TRUE :4255      TRUE :377
#  Layer2-Layer4   Layer2-Layer5   Layer2-Layer6   Layer3-Layer4
#  Mode :logical   Mode :logical   Mode :logical   Mode :logical
#  FALSE:20047     FALSE:20026     FALSE:19884     FALSE:22004
#  TRUE :2284      TRUE :2305      TRUE :2447      TRUE :327
#  Layer3-Layer5   Layer3-Layer6   Layer4-Layer5   Layer4-Layer6
#  Mode :logical   Mode :logical   Mode :logical   Mode :logical
#  FALSE:21389     FALSE:20579     FALSE:22039     FALSE:20586
#  TRUE :942       TRUE :1752      TRUE :292       TRUE :1745
#  Layer5-Layer6
#  Mode :logical
#  FALSE:21816
#  TRUE :515
sort(colSums(fdrs_contrasts < 0.05))
# Layer4-Layer5 Layer3-Layer4 Layer2-Layer3 Layer5-Layer6 Layer3-Layer5
#           292           327           377           515           942
# Layer4-Layer6 Layer3-Layer6 Layer2-Layer4 Layer2-Layer5 Layer2-Layer6
#          1745          1752          2284          2305          2447
# Layer1-Layer3 Layer1-Layer2 Layer1-Layer6 Layer1-Layer5 Layer1-Layer4
#          3566          3686          4255          4638          4677
#     WM-Layer1     WM-Layer6     WM-Layer5     WM-Layer2     WM-Layer4
#          5667          6686          7979          7992          8458
#     WM-Layer3
#          8500


## Make some supplementary tables
f_merge <- function(p, fdr, t) {
    colnames(p) <- paste0('p_value_', colnames(p))
    colnames(fdr) <- paste0('fdr_', colnames(fdr))
    colnames(t) <- paste0('t_stat_', colnames(t))
    res <- as.data.frame(cbind(t, p, fdr))
    res$ensembl <- rownames(res)
    ## Check it's all in order
    stopifnot(identical(rownames(res), rownames(sce_layer)))
    res$gene <- rowData(sce_layer)$gene_name
    rownames(res) <- NULL
    return(res)
}

results_specificity <-
    f_merge(p = pvals0_contrasts, fdr = fdrs0_contrasts, t = t0_contrasts)
head(results_specificity)
#    t_stat_WM t_stat_Layer1 t_stat_Layer2 t_stat_Layer3 t_stat_Layer4
# 1 -0.6344143    -1.0321320    0.17815008   -0.72835965    1.56703859
# 2 -2.4758891     1.2232062   -0.87337451    1.93793650    1.33150141
# 3 -3.0079360    -0.8564572    2.13358520    0.48741121    0.35212807
# 4 -1.2916584    -0.9494234   -0.94854397    0.56378302   -0.11206713
# 5  2.3175897     0.6156900    0.11274780   -0.09907566   -0.03376771
# 6 -2.2686017    -0.6536163   -0.08615251    1.84786166    0.77710957
#   t_stat_Layer5 t_stat_Layer6  p_value_WM p_value_Layer1 p_value_Layer2
# 1    -0.2202707     0.7438713 0.527700348      0.3052551     0.85907467
# 2     0.4773214    -1.6152865 0.015497447      0.2249981     0.38518446
# 3     0.6071363     0.3832779 0.003557367      0.3944138     0.03607273
# 4     1.1536114     1.2634726 0.200358650      0.3453884     0.34583092
# 5    -0.3434730    -2.4629827 0.023143222      0.5399225     0.91052499
# 6     0.0355838     0.1945122 0.026108214      0.5153146     0.93156963
#   p_value_Layer3 p_value_Layer4 p_value_Layer5 p_value_Layer6     fdr_WM
# 1     0.46861049      0.1212213      0.8262455     0.45922628 0.63711486
# 2     0.05630996      0.1869669      0.6344905     0.11035355 0.03959651
# 3     0.62735659      0.7257077      0.5455535     0.70257361 0.01107944
# 4     0.57454588      0.9110629      0.2522414     0.21024494 0.31550730
# 5     0.92133658      0.9731501      0.7321823     0.01601949 0.05551142
# 6     0.06847596      0.4394838      0.9717067     0.84628893 0.06148729
#   fdr_Layer1 fdr_Layer2 fdr_Layer3 fdr_Layer4 fdr_Layer5 fdr_Layer6
# 1  0.5644497  0.9418694  0.8284720  0.4139767  0.9596911  0.8115244
# 2  0.4828399  0.6944277  0.3831776  0.4938944  0.9051814  0.5481106
# 3  0.6380674  0.2127681  0.9022535  0.8830293  0.8698217  0.9150742
# 4  0.6001298  0.6683198  0.8845131  0.9674045  0.7017599  0.6698502
# 5  0.7356770  0.9635319  0.9848034  0.9910492  0.9343088  0.2431891
# 6  0.7192630  0.9733708  0.4139515  0.7028656  0.9936069  0.9633326
#           ensembl        gene
# 1 ENSG00000243485 MIR1302-2HG
# 2 ENSG00000238009  AL627309.1
# 3 ENSG00000237491  AL669831.5
# 4 ENSG00000177757      FAM87B
# 5 ENSG00000225880   LINC00115
# 6 ENSG00000230368      FAM41C


results_pairwise <-
    f_merge(p = pvals_contrasts, fdr = fdrs_contrasts, t = eb_contrasts$t)
colnames(results_pairwise)
#  [1] "t_stat_WM-Layer1"      "t_stat_WM-Layer2"      "t_stat_WM-Layer3"
#  [4] "t_stat_WM-Layer4"      "t_stat_WM-Layer5"      "t_stat_WM-Layer6"
#  [7] "t_stat_Layer1-Layer2"  "t_stat_Layer1-Layer3"  "t_stat_Layer1-Layer4"
# [10] "t_stat_Layer1-Layer5"  "t_stat_Layer1-Layer6"  "t_stat_Layer2-Layer3"
# [13] "t_stat_Layer2-Layer4"  "t_stat_Layer2-Layer5"  "t_stat_Layer2-Layer6"
# [16] "t_stat_Layer3-Layer4"  "t_stat_Layer3-Layer5"  "t_stat_Layer3-Layer6"
# [19] "t_stat_Layer4-Layer5"  "t_stat_Layer4-Layer6"  "t_stat_Layer5-Layer6"
# [22] "p_value_WM-Layer1"     "p_value_WM-Layer2"     "p_value_WM-Layer3"
# [25] "p_value_WM-Layer4"     "p_value_WM-Layer5"     "p_value_WM-Layer6"
# [28] "p_value_Layer1-Layer2" "p_value_Layer1-Layer3" "p_value_Layer1-Layer4"
# [31] "p_value_Layer1-Layer5" "p_value_Layer1-Layer6" "p_value_Layer2-Layer3"
# [34] "p_value_Layer2-Layer4" "p_value_Layer2-Layer5" "p_value_Layer2-Layer6"
# [37] "p_value_Layer3-Layer4" "p_value_Layer3-Layer5" "p_value_Layer3-Layer6"
# [40] "p_value_Layer4-Layer5" "p_value_Layer4-Layer6" "p_value_Layer5-Layer6"
# [43] "fdr_WM-Layer1"         "fdr_WM-Layer2"         "fdr_WM-Layer3"
# [46] "fdr_WM-Layer4"         "fdr_WM-Layer5"         "fdr_WM-Layer6"
# [49] "fdr_Layer1-Layer2"     "fdr_Layer1-Layer3"     "fdr_Layer1-Layer4"
# [52] "fdr_Layer1-Layer5"     "fdr_Layer1-Layer6"     "fdr_Layer2-Layer3"
# [55] "fdr_Layer2-Layer4"     "fdr_Layer2-Layer5"     "fdr_Layer2-Layer6"
# [58] "fdr_Layer3-Layer4"     "fdr_Layer3-Layer5"     "fdr_Layer3-Layer6"
# [61] "fdr_Layer4-Layer5"     "fdr_Layer4-Layer6"     "fdr_Layer5-Layer6"
# [64] "ensembl"               "gene"


## Match the colnames to the new style
f_rename <- function(x, old, new = old) {
    old_patt <- paste0('_', old, '$')
    i <- grep(old_patt, colnames(x))
    tmp <- gsub(old_patt, '', colnames(x)[i])
    tmp <- paste0(new, '_', tmp)
    colnames(x)[i] <- tmp
    return(x)
}
results_anova <-
    f_rename(f_rename(f_rename(
        f_rename(f_stats, 'f', 'f_stat'), 'p_value'
    ), 'fdr'), 'Amean')
head(results_anova)
#   f_stat_full p_value_full     fdr_full Amean_full f_stat_noWM p_value_noWM
# 1    1.126228 3.565528e-01 3.679551e-01  0.1422407    1.094988 3.758663e-01
# 2    8.282190 2.294975e-07 3.118859e-07  0.8383094    8.080743 2.028620e-06
# 3   90.122437 7.145239e-33 1.319553e-32  3.5454133  222.281074 3.587143e-39
# 4    1.447935 2.000407e-01 2.120027e-01  0.1048252    1.314147 2.648163e-01
# 5   10.377382 6.555601e-09 9.162741e-09  1.2135801   10.136817 9.944087e-08
# 6   10.650917 4.217147e-09 5.914653e-09  1.3447335   12.470669 4.417884e-09
#       fdr_noWM Amean_noWM         ensembl        gene
# 1 3.875101e-01  0.1565079 ENSG00000243485 MIR1302-2HG
# 2 2.732439e-06  0.9362990 ENSG00000238009  AL627309.1
# 3 7.312807e-39  3.6730061 ENSG00000237491  AL669831.5
# 4 2.778953e-01  0.1244799 ENSG00000177757      FAM87B
# 5 1.372529e-07  1.1121321 ENSG00000225880   LINC00115
# 6 6.250761e-09  1.4543898 ENSG00000230368      FAM41C

## Save for later
save(results_anova, results_specificity, results_pairwise, file = 'rda/modeling_results.Rdata')


## Make a ggplot-version of
## https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/Layer_Guesses/layer_specificity_fstats.R#L93-L125
## with colors by some variables

anova_df <- results_anova
## Compute the percent of WM expression
ix_wm <- which(sce_layer$layer_guess == 'WM')
anova_df$expr_total <- rowSums(assays(sce_layer)$counts)
anova_df$expr_WM <- rowSums(assays(sce_layer)$counts[, ix_wm])
anova_df$expr_WM_ratio <- anova_df$expr_WM / anova_df$expr_total


pdf(
    'pdf/layer_specificity_full_vs_noWM_WMratio.pdf',
    useDingbats = FALSE,
    width = 10
)
ggplot(anova_df,
    aes(x = f_stat_full, y = f_stat_noWM, color = expr_WM_ratio)) +
    geom_point() + theme_bw(base_size = 20) +
    scale_color_gradientn(name = 'WM ratio', colors = viridis(21)) +
    geom_abline(
        intercept = 0,
        slope = 1,
        colour = 'red',
        linetype = 2
    ) +
    xlab('F-stats (WM + L1 through L6)') +
    ylab('F-stats (L1 through L6 only)')
ggplot(anova_df,
    aes(
        x = -log10(p_value_full),
        y = -log10(p_value_noWM),
        color = expr_WM_ratio
    )) +
    geom_point() + theme_bw(base_size = 20) +
    scale_color_gradientn(name = 'WM ratio', colors = viridis(21)) +
    geom_abline(
        intercept = 0,
        slope = 1,
        colour = 'red',
        linetype = 2
    )  +
    xlab('-log10 p-value (WM + L1 through L6)') +
    ylab('-log10 p-value (L1 through L6 only)')
dev.off()

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
#  date     2020-02-13
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  assertthat             0.2.1     2019-03-21 [2] CRAN (R 3.6.1)
#  backports              1.1.5     2019-10-02 [1] CRAN (R 3.6.1)
#  Biobase              * 2.46.0    2019-10-29 [2] Bioconductor
#  BiocGenerics         * 0.32.0    2019-10-29 [1] Bioconductor
#  BiocParallel         * 1.20.1    2019-12-21 [1] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.6.1)
#  cellranger             1.1.0     2016-07-27 [1] CRAN (R 3.6.1)
#  cli                    2.0.0     2019-12-09 [1] CRAN (R 3.6.1)
#  colorout             * 1.2-2     2019-10-31 [1] Github (jalvesaq/colorout@641ed38)
#  colorspace             1.4-1     2019-03-18 [2] CRAN (R 3.6.1)
#  cowplot              * 1.0.0     2019-07-11 [1] CRAN (R 3.6.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.6.1)
#  DelayedArray         * 0.12.0    2019-10-29 [2] Bioconductor
#  digest                 0.6.23    2019-11-23 [1] CRAN (R 3.6.1)
#  dplyr                  0.8.3     2019-07-04 [1] CRAN (R 3.6.1)
#  fansi                  0.4.0     2018-10-05 [1] CRAN (R 3.6.1)
#  GenomeInfoDb         * 1.22.0    2019-10-29 [1] Bioconductor
#  GenomeInfoDbData       1.2.2     2019-10-28 [2] Bioconductor
#  GenomicRanges        * 1.38.0    2019-10-29 [1] Bioconductor
#  ggplot2              * 3.2.1     2019-08-10 [1] CRAN (R 3.6.1)
#  glue                   1.3.1     2019-03-12 [1] CRAN (R 3.6.1)
#  gridExtra            * 2.3       2017-09-09 [2] CRAN (R 3.6.1)
#  gtable                 0.3.0     2019-03-25 [2] CRAN (R 3.6.1)
#  here                 * 0.1       2017-05-28 [1] CRAN (R 3.6.1)
#  htmltools              0.4.0     2019-10-04 [1] CRAN (R 3.6.1)
#  htmlwidgets            1.5.1     2019-10-08 [1] CRAN (R 3.6.1)
#  httpuv                 1.5.2     2019-09-11 [1] CRAN (R 3.6.1)
#  IRanges              * 2.20.1    2019-11-20 [1] Bioconductor
#  jsonlite               1.6       2018-12-07 [2] CRAN (R 3.6.1)
#  labeling               0.3       2014-08-23 [2] CRAN (R 3.6.1)
#  later                  1.0.0     2019-10-04 [1] CRAN (R 3.6.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.6.1)
#  lazyeval               0.2.2     2019-03-15 [2] CRAN (R 3.6.1)
#  limma                * 3.42.0    2019-10-29 [1] Bioconductor
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.6.1)
#  Matrix                 1.2-17    2019-03-22 [3] CRAN (R 3.6.1)
#  matrixStats          * 0.55.0    2019-09-07 [1] CRAN (R 3.6.1)
#  mime                   0.8       2019-12-19 [1] CRAN (R 3.6.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.6.1)
#  pillar                 1.4.3     2019-12-20 [1] CRAN (R 3.6.1)
#  pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 3.6.1)
#  plyr                   1.8.4     2016-06-08 [2] CRAN (R 3.6.1)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.6.1)
#  Polychrome           * 1.2.3     2019-08-01 [1] CRAN (R 3.6.1)
#  promises               1.1.0     2019-10-04 [1] CRAN (R 3.6.1)
#  purrr                  0.3.3     2019-10-18 [2] CRAN (R 3.6.1)
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.6.1)
#  RColorBrewer           1.1-2     2014-12-07 [2] CRAN (R 3.6.1)
#  Rcpp                   1.0.3     2019-11-08 [1] CRAN (R 3.6.1)
#  RCurl                  1.95-4.12 2019-03-04 [2] CRAN (R 3.6.1)
#  readxl               * 1.3.1     2019-03-13 [2] CRAN (R 3.6.1)
#  reshape2               1.4.3     2017-12-11 [2] CRAN (R 3.6.1)
#  rlang                  0.4.2     2019-11-23 [1] CRAN (R 3.6.1)
#  rmote                * 0.3.4     2019-10-31 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot              1.3-2     2018-01-03 [2] CRAN (R 3.6.1)
#  S4Vectors            * 0.24.1    2019-12-01 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.6.1)
#  scatterplot3d          0.3-41    2018-03-14 [1] CRAN (R 3.6.1)
#  servr                  0.15      2019-08-07 [1] CRAN (R 3.6.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.6.1)
#  SingleCellExperiment * 1.8.0     2019-10-29 [2] Bioconductor
#  stringi                1.4.3     2019-03-12 [2] CRAN (R 3.6.1)
#  stringr                1.4.0     2019-02-10 [1] CRAN (R 3.6.1)
#  SummarizedExperiment * 1.16.1    2019-12-19 [1] Bioconductor
#  tibble                 2.1.3     2019-06-06 [1] CRAN (R 3.6.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.6.1)
#  viridisLite          * 0.3.0     2018-02-01 [2] CRAN (R 3.6.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.6.1)
#  xfun                   0.11      2019-11-12 [1] CRAN (R 3.6.1)
#  XVector                0.26.0    2019-10-29 [1] Bioconductor
#  zlibbioc               1.32.0    2019-10-29 [2] Bioconductor
#
# [1] /users/lcollado/R/3.6.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/library

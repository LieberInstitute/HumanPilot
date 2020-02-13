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
#  gtable                 0.3.0     2019-03-25 [2] CRAN (R 3.6.1)
#  here                 * 0.1       2017-05-28 [1] CRAN (R 3.6.1)
#  htmltools              0.4.0     2019-10-04 [1] CRAN (R 3.6.1)
#  htmlwidgets            1.5.1     2019-10-08 [1] CRAN (R 3.6.1)
#  httpuv                 1.5.2     2019-09-11 [1] CRAN (R 3.6.1)
#  IRanges              * 2.20.1    2019-11-20 [1] Bioconductor
#  jsonlite               1.6       2018-12-07 [2] CRAN (R 3.6.1)
#  later                  1.0.0     2019-10-04 [1] CRAN (R 3.6.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.6.1)
#  lazyeval               0.2.2     2019-03-15 [2] CRAN (R 3.6.1)
#  limma                * 3.42.0    2019-10-29 [1] Bioconductor
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.6.1)
#  Matrix                 1.2-17    2019-03-22 [3] CRAN (R 3.6.1)
#  matrixStats          * 0.55.0    2019-09-07 [1] CRAN (R 3.6.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.6.1)
#  pillar                 1.4.3     2019-12-20 [1] CRAN (R 3.6.1)
#  pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 3.6.1)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.6.1)
#  Polychrome           * 1.2.3     2019-08-01 [1] CRAN (R 3.6.1)
#  promises               1.1.0     2019-10-04 [1] CRAN (R 3.6.1)
#  purrr                  0.3.3     2019-10-18 [2] CRAN (R 3.6.1)
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.6.1)
#  RColorBrewer           1.1-2     2014-12-07 [2] CRAN (R 3.6.1)
#  Rcpp                   1.0.3     2019-11-08 [1] CRAN (R 3.6.1)
#  RCurl                  1.95-4.12 2019-03-04 [2] CRAN (R 3.6.1)
#  readxl               * 1.3.1     2019-03-13 [2] CRAN (R 3.6.1)
#  rlang                  0.4.2     2019-11-23 [1] CRAN (R 3.6.1)
#  rmote                * 0.3.4     2019-10-31 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot              1.3-2     2018-01-03 [2] CRAN (R 3.6.1)
#  S4Vectors            * 0.24.1    2019-12-01 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.6.1)
#  scatterplot3d          0.3-41    2018-03-14 [1] CRAN (R 3.6.1)
#  servr                  0.15      2019-08-07 [1] CRAN (R 3.6.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.6.1)
#  SingleCellExperiment * 1.8.0     2019-10-29 [2] Bioconductor
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

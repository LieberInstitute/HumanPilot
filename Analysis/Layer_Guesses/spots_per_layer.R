library('SingleCellExperiment')
library('here')
library('dplyr')
library('sessioninfo')

## Load data
load(here(
    'Analysis',
    'Human_DLPFC_Visium_processedData_sce_scran.Rdata'
))


tab <-
    read.csv(
        dir(
            here('Analysis', 'Layer_Guesses', 'Second_Round'),
            pattern = 'Combined2',
            full.names = TRUE
        ),
        header = TRUE,
        stringsAsFactors = FALSE,
        na.strings = ''
    )

tab$key <- paste0(tab$sample_name,
    '_',
    tab$spot_name)


## From below, the guesses are for samples 3,4 + 7,8 + 11,12
## (the ones 300 microns away)
unique(sce$sample_name)
# [1] 151507 151508 151509 151510 151669 151670 151671 151672 151673 151674 151675 151676
unique(tab$sample_name)
# [1] 151509 151510 151671 151672 151675 151676

## Spots per layer
with(tab, tapply(layer, sample_name, table))

## tidyverse-way
spots_layer <-
    group_by(tab, sample_name, layer) %>% summarize(spots_per_layer = n())

## Add the subject
spots_layer$subject <-
    sce$subject[match(spots_layer$sample_name, sce$sample_name)]

## Save for later
save(spots_layer,
    file = here('Analysis', 'Layer_Guesses', 'spots_layer.Rdata'))

## Global mean: 641.1 spots per layer across all 6 slides
## and across all layers
summary(spots_layer$spots_per_layer)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 166.0   292.8   533.0   641.1   753.0  1918.0

## Mean number of spots per layer for each of the 6 slides
summarize(spots_layer, mean_spots = mean(spots_per_layer))
# # A tibble: 6 x 2
# sample_name mean_spots
# <int>      <dbl>
#     1      151509       684
# 2      151510       656.
# 3      151671       819.
# 4      151672       778.
# 5      151675       509.
# 6      151676       490.

## Mean number of spots per layer for each of the different layers
## across all 6 images
group_by(spots_layer, layer) %>% summarize(mean_spots = mean(spots_per_layer))
# # A tibble: 8 x 2
# layer     mean_spots
# <chr>          <dbl>
#     1 Layer 1         746.
# 2 Layer 2         445.
# 3 Layer 2/3      1746.
# 4 Layer 3        1316.
# 5 Layer 4         294.
# 6 Layer 5         584.
# 7 Layer 6         531.
# 8 WM              397.

## Mean number of spots per layer across each of the 3 subjects
## (so data from 2 slides for each subject)
group_by(spots_layer, subject, layer) %>% summarize(mean_spots = mean(spots_per_layer))
# # A tibble: 19 x 3
# # Groups:   subject [3]
# subject layer     mean_spots
# <chr>   <chr>          <dbl>
#     1 Br5292  Layer 1        1184.
# 2 Br5292  Layer 2         626
# 3 Br5292  Layer 3        1829
# 4 Br5292  Layer 4         344.
# 5 Br5292  Layer 5         336.
# 6 Br5292  Layer 6         197
# 7 Br5292  WM              175
# 8 Br5595  Layer 2/3      1746.
# 9 Br5595  Layer 4         274.
# 10 Br5595  Layer 5         724.
# 11 Br5595  Layer 6         821
# 12 Br5595  WM              424
# 13 Br8100  Layer 1         308.
# 14 Br8100  Layer 2         264.
# 15 Br8100  Layer 3         804.
# 16 Br8100  Layer 4         264.
# 17 Br8100  Layer 5         690.
# 18 Br8100  Layer 6         574.
# 19 Br8100  WM              592.

## Also check the number of cells per spot

## First globally
summary(sce$cell_count)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.000   1.000   3.000   3.281   4.000  27.000

## Then by sample
group_by(as.data.frame(colData(sce)), sample_name) %>% summarize(mean_cells = mean(cell_count))
# # A tibble: 12 x 2
# sample_name mean_cells
# <int>      <dbl>
#     1      151507       2.20
# 2      151508       3.31
# 3      151509       2.99
# 4      151510       3.20
# 5      151669       2.63
# 6      151670       5.59
# 7      151671       2.73
# 8      151672       2.07
# 9      151673       4.53
# 10      151674       3.97
# 11      151675       3.46
# 12      151676       3.24

## Next by subject
group_by(as.data.frame(colData(sce)), subject) %>% summarize(mean_cells = mean(cell_count))
# # A tibble: 3 x 2
# subject mean_cells
# <chr>        <dbl>
#     1 Br5292        2.94
# 2 Br5595        3.19
# 3 Br8100        3.81

## Finally by subject and position
group_by(as.data.frame(colData(sce)), subject, position) %>% summarize(mean_cells = mean(cell_count))
# # A tibble: 6 x 3
# # Groups:   subject [3]
# subject position mean_cells
# <chr>   <chr>         <dbl>
#     1 Br5292  0              2.77
# 2 Br5292  300            3.10
# 3 Br5595  0              4.08
# 4 Br5595  300            2.40
# 5 Br8100  0              4.25
# 6 Br8100  300            3.35

## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 3.6.2 (2019-12-12)
# os       Windows 10 x64
# system   x86_64, mingw32
# ui       RStudio
# language (EN)
# collate  English_United States.1252
# ctype    English_United States.1252
# tz       America/Mexico_City
# date     2020-01-07
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version   date       lib source
# assertthat             0.2.1     2019-03-21 [1] CRAN (R 3.6.1)
# backports              1.1.5     2019-10-02 [1] CRAN (R 3.6.1)
# Biobase              * 2.46.0    2019-10-29 [1] Bioconductor
# BiocGenerics         * 0.32.0    2019-10-29 [1] Bioconductor
# BiocParallel         * 1.20.1    2019-12-21 [1] Bioconductor
# bitops                 1.0-6     2013-08-17 [1] CRAN (R 3.6.0)
# cli                    2.0.0     2019-12-09 [1] CRAN (R 3.6.1)
# crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.6.1)
# DelayedArray         * 0.12.1    2019-12-17 [1] Bioconductor
# dplyr                * 0.8.3     2019-07-04 [1] CRAN (R 3.6.1)
# fansi                  0.4.0     2018-10-05 [1] CRAN (R 3.6.1)
# GenomeInfoDb         * 1.22.0    2019-10-29 [1] Bioconductor
# GenomeInfoDbData       1.2.2     2019-12-18 [1] Bioconductor
# GenomicRanges        * 1.38.0    2019-10-29 [1] Bioconductor
# glue                   1.3.1     2019-03-12 [1] CRAN (R 3.6.1)
# here                 * 0.1       2017-05-28 [1] CRAN (R 3.6.1)
# IRanges              * 2.20.1    2019-11-20 [1] Bioconductor
# lattice                0.20-38   2018-11-04 [2] CRAN (R 3.6.2)
# magrittr               1.5       2014-11-22 [1] CRAN (R 3.6.1)
# Matrix                 1.2-18    2019-11-27 [2] CRAN (R 3.6.2)
# matrixStats          * 0.55.0    2019-09-07 [1] CRAN (R 3.6.1)
# packrat                0.5.0     2018-11-14 [1] CRAN (R 3.6.1)
# pillar                 1.4.3     2019-12-20 [1] CRAN (R 3.6.2)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 3.6.1)
# purrr                  0.3.3     2019-10-18 [1] CRAN (R 3.6.1)
# R6                     2.4.1     2019-11-12 [1] CRAN (R 3.6.1)
# Rcpp                   1.0.3     2019-11-08 [1] CRAN (R 3.6.1)
# RCurl                  1.95-4.12 2019-03-04 [1] CRAN (R 3.6.0)
# rlang                  0.4.2     2019-11-23 [1] CRAN (R 3.6.1)
# rprojroot              1.3-2     2018-01-03 [1] CRAN (R 3.6.1)
# rstudioapi             0.10      2019-03-19 [1] CRAN (R 3.6.1)
# S4Vectors            * 0.24.1    2019-12-01 [1] Bioconductor
# sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.6.1)
# SingleCellExperiment * 1.8.0     2019-10-29 [1] Bioconductor
# SummarizedExperiment * 1.16.1    2019-12-20 [1] Bioconductor
# tibble                 2.1.3     2019-06-06 [1] CRAN (R 3.6.1)
# tidyselect             0.2.5     2018-10-11 [1] CRAN (R 3.6.1)
# utf8                   1.1.4     2018-05-24 [1] CRAN (R 3.6.1)
# vctrs                  0.2.1     2019-12-17 [1] CRAN (R 3.6.2)
# withr                  2.1.2     2018-03-15 [1] CRAN (R 3.6.1)
# XVector                0.26.0    2019-10-29 [1] Bioconductor
# zeallot                0.1.0     2018-01-28 [1] CRAN (R 3.6.1)
# zlibbioc               1.32.0    2019-10-29 [1] Bioconductor
#
# [1] D:/Documents/R/win-library/3.6
# [2] D:/R/R-3.6.2/library

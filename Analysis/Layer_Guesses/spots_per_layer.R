library('SingleCellExperiment')
library('here')
library('dplyr')
library('sessioninfo')

## Load data
load(here(
    'Analysis',
    'Human_DLPFC_Visium_processedData_sce_scran.Rdata'
))

## For plotting
source(here('Analysis', 'spatialLIBD_global_plot_code.R'))
genes <- paste0(rowData(sce)$gene_name, '; ', rowData(sce)$gene_id)

## Merge the first and second round merged guesses
tab <-
    read.csv(
        dir(
            here('Analysis', 'Layer_Guesses', 'First_Round'),
            pattern = 'Merged',
            full.names = TRUE
        ),
        header = TRUE,
        stringsAsFactors = FALSE,
        na.strings = ''
    )

tab$key <- paste0(tab$sample_name,
    '_',
    tab$spot_name)

tab2 <-
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

tab2$key <- paste0(tab2$sample_name,
    '_',
    tab2$spot_name)

stopifnot(!any(tab$key %in% tab2$key))
stopifnot(!any(tab2$key %in% tab$key))

## Save for other work
layer_guess_tab <- rbind(tab, tab2)
rownames(layer_guess_tab) <- NULL
save(layer_guess_tab,
    file = here('Analysis', 'Layer_Guesses', 'rda', 'layer_guess_tab.Rdata'))

## Find which ones are missing:
m <- match(sce$key, layer_guess_tab$key)
options(width = 100)
addmargins(table(sce$sample_name[is.na(m)]))
# 151507 151508 151509 151510 151669 151670 151671 151672 151673 151674 151675 151676    Sum
#      5      3      1     39     25     14     17    127     28     38     26     29    352

## Spots per layer
with(layer_guess_tab, tapply(layer, sample_name, table))

## tidyverse-way
spots_layer <-
    group_by(tab, sample_name, layer) %>% summarize(spots_per_layer = n())

## Add the subject
spots_layer$subject <-
    sce$subject[match(spots_layer$sample_name, sce$sample_name)]

## Save for later
save(spots_layer,
    file = here('Analysis', 'Layer_Guesses', 'rda', 'spots_layer.Rdata'))

## Global mean: 604.4 spots per layer across all slides
## and across all layers
summary(spots_layer$spots_per_layer)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 200.0   297.5   498.0   604.4   687.8  2175.0

## Mean number of spots per layer for each of the 12 slides
summarize(spots_layer, mean_spots = mean(spots_per_layer))
# # A tibble: 6 x 2
#   sample_name mean_spots
#         <int>      <dbl>
# 1      151507       603
# 2      151508       626.
# 3      151669       727.
# 4      151670       697.
# 5      151673       516.
# 6      151674       519.

## Mean number of spots per layer for each of the different layers
## across all 12 images
group_by(spots_layer, layer) %>% summarize(mean_spots = mean(spots_per_layer))
# # A tibble: 8 x 2
#   layer     mean_spots
#   <chr>          <dbl>
# 1 Layer 1         584
# 2 Layer 2         269.
# 3 Layer 2/3      2158
# 4 Layer 3        1128.
# 5 Layer 4         297
# 6 Layer 5         633.
# 7 Layer 6         503.
# 8 WM              355.

## Mean number of spots per layer across each of the 3 subjects
group_by(spots_layer, subject, layer) %>% summarize(mean_spots = mean(spots_per_layer))
# # A tibble: 19 x 3
# # Groups:   subject [3]
#    subject layer     mean_spots
#    <chr>   <chr>          <dbl>
#  1 Br5292  Layer 1         842.
#  2 Br5292  Layer 2         300
#  3 Br5292  Layer 3        1300
#  4 Br5292  Layer 4         371
#  5 Br5292  Layer 5         706
#  6 Br5292  Layer 6         506.
#  7 Br5292  WM              277
#  8 Br5595  Layer 2/3      2158
#  9 Br5595  Layer 4         288.
# 10 Br5595  Layer 5         546.
# 11 Br5595  Layer 6         350.
# 12 Br5595  WM              220.
# 13 Br8100  Layer 1         326.
# 14 Br8100  Layer 2         238.
# 15 Br8100  Layer 3         956.
# 16 Br8100  Layer 4         232.
# 17 Br8100  Layer 5         647
# 18 Br8100  Layer 6         653
# 19 Br8100  WM              569

## Also check the number of cells per spot

## First globally
summary(sce$cell_count)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.000   1.000   3.000   3.281   4.000  27.000

sce_image_grid_gene(
    sce,
    geneid = 'cell_count',
    spatial = TRUE,
    minCount = 0,
    pdf_file = 'pdf/spot_cell_count.pdf'
)


## Then by sample
group_by(as.data.frame(colData(sce)), sample_name) %>% summarize(mean_cells = mean(cell_count))
# # A tibble: 12 x 2
#    sample_name mean_cells
#          <int>      <dbl>
#  1      151507       2.20
#  2      151508       3.31
#  3      151509       2.99
#  4      151510       3.20
#  5      151669       2.63
#  6      151670       5.59
#  7      151671       2.73
#  8      151672       2.07
#  9      151673       4.53
# 10      151674       3.97
# 11      151675       3.46
# 12      151676       3.24

group_by(as.data.frame(colData(sce)), sample_name) %>% summarize(mean_cells = mean(cell_count)) %>% summary()
# sample_name       mean_cells
# Min.   :151507   Min.   :2.066
# 1st Qu.:151510   1st Qu.:2.708
# Median :151670   Median :3.222
# Mean   :151618   Mean   :3.328
# 3rd Qu.:151673   3rd Qu.:3.586
# Max.   :151676   Max.   :5.595

layer_guess_tab$cell_count = sce$cell_count[match(layer_guess_tab$key, sce$key)]

layer_guess_tab$layer_guess = layer_guess_tab$layer
layer_guess_tab$layer_guess[layer_guess_tab$layer_guess == 'Layer 2/3'] <- 'Layer 3'

## This is supp table 2
group_by(layer_guess_tab, layer_guess) %>%
    summarize(prop0 = mean(cell_count == 0),
        prop1 = mean(cell_count == 1))
# # A tibble: 8 x 3
  # layer_guess  prop0  prop1
  # <chr>        <dbl>  <dbl>
# 1 Layer 1     0.334  0.217
# 2 Layer 2     0.142  0.179
# 3 Layer 3     0.104  0.176
# 4 Layer 4     0.0313 0.118
# 5 Layer 5     0.0344 0.119
# 6 Layer 6     0.0484 0.128
# 7 WM          0.0385 0.0594


group_by(layer_guess_tab, layer) %>%
    summarize(prop0 = mean(cell_count == 0),
        prop1 = mean(cell_count == 1)) %>% summary()
#    layer               prop0             prop1
# Length:8           Min.   :0.03129   Min.   :0.05937
# Class :character   1st Qu.:0.03751   1st Qu.:0.11854
# Mode  :character   Median :0.07144   Median :0.14696
#                    Mean   :0.10442   Mean   :0.14678
#                    3rd Qu.:0.11949   3rd Qu.:0.18122
#                    Max.   :0.33446   Max.   :0.21702

group_by(layer_guess_tab, sample_name) %>%
    summarize(prop0 = mean(cell_count == 0),
        prop1 = mean(cell_count == 1))
# # A tibble: 12 x 3
#    sample_name  prop0  prop1
#          <int>  <dbl>  <dbl>
#  1      151507 0.243  0.201
#  2      151508 0.132  0.139
#  3      151509 0.216  0.132
#  4      151510 0.101  0.133
#  5      151669 0.0589 0.176
#  6      151670 0.0227 0.0672
#  7      151671 0.0584 0.164
#  8      151672 0.120  0.263
#  9      151673 0.0233 0.0584
# 10      151674 0.0413 0.114
# 11      151675 0.0645 0.160
# 12      151676 0.0851 0.187

group_by(layer_guess_tab, sample_name) %>%
    summarize(prop0 = mean(cell_count == 0),
        prop1 = mean(cell_count == 1)) %>% summary()
#  sample_name         prop0             prop1
# Min.   :151507   Min.   :0.02268   Min.   :0.05843
# 1st Qu.:151510   1st Qu.:0.05411   1st Qu.:0.12772
# Median :151670   Median :0.07480   Median :0.14957
# Mean   :151618   Mean   :0.09714   Mean   :0.14956
# 3rd Qu.:151673   3rd Qu.:0.12301   3rd Qu.:0.17885
# Max.   :151676   Max.   :0.24331   Max.   :0.26312

group_by(layer_guess_tab, layer, sample_name) %>%
    summarize(prop0 = mean(cell_count == 0),
        prop1 = mean(cell_count == 1))
# # A tibble: 76 x 4
# # Groups:   layer [8]
#    layer   sample_name prop0 prop1
#    <chr>         <int> <dbl> <dbl>
#  1 Layer 1      151507 0.704 0.165
#  2 Layer 1      151508 0.486 0.260
#  3 Layer 1      151509 0.168 0.197
#  4 Layer 1      151510 0.312 0.261
#  5 Layer 1      151673 0.143 0.172
#  6 Layer 1      151674 0.166 0.192
#  7 Layer 1      151675 0.134 0.186
#  8 Layer 1      151676 0.242 0.249
#  9 Layer 2      151507 0.364 0.305
# 10 Layer 2      151508 0.163 0.224
# # … with 66 more rows

group_by(layer_guess_tab, layer, sample_name) %>%
    summarize(prop0 = mean(cell_count == 0),
        prop1 = mean(cell_count == 1)) %>% summary()
#    layer            sample_name         prop0             prop1
# Length:76          Min.   :151507   Min.   :0.00000   Min.   :0.0064
# Class :character   1st Qu.:151509   1st Qu.:0.01287   1st Qu.:0.0524
# Mode  :character   Median :151670   Median :0.04030   Median :0.1284
#                    Mean   :151612   Mean   :0.09198   Mean   :0.1342
#                    3rd Qu.:151674   3rd Qu.:0.10067   3rd Qu.:0.1955
#                    Max.   :151676   Max.   :0.74699   Max.   :0.3049

## Next by subject
group_by(as.data.frame(colData(sce)), subject) %>% summarize(mean_cells = mean(cell_count))
# # A tibble: 3 x 2
#   subject mean_cells
#   <chr>        <dbl>
# 1 Br5292        2.94
# 2 Br5595        3.19
# 3 Br8100        3.81

group_by(as.data.frame(colData(sce)), subject) %>% summarize(mean_cells = mean(cell_count)) %>%  summary()
#   subject            mean_cells
# Length:3           Min.   :2.939
# Class :character   1st Qu.:3.064
# Mode  :character   Median :3.188
#                    Mean   :3.312
#                    3rd Qu.:3.499
#                    Max.   :3.809

## Finally by subject and position
group_by(as.data.frame(colData(sce)), subject, position) %>% summarize(mean_cells = mean(cell_count))
# # A tibble: 6 x 3
# # Groups:   subject [3]
#   subject position mean_cells
#   <chr>   <chr>         <dbl>
# 1 Br5292  0              2.77
# 2 Br5292  300            3.10
# 3 Br5595  0              4.08
# 4 Br5595  300            2.40
# 5 Br8100  0              4.25
# 6 Br8100  300            3.35

group_by(as.data.frame(colData(sce)), subject, position) %>% summarize(mean_cells = mean(cell_count)) %>% summary()
#   subject            position           mean_cells
# Length:6           Length:6           Min.   :2.403
# Class :character   Class :character   1st Qu.:2.849
# Mode  :character   Mode  :character   Median :3.224
#                                       Mean   :3.324
#                                       3rd Qu.:3.898
#                                       Max.   :4.250

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
#  date     2020-01-16
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  assertthat             0.2.1     2019-03-21 [2] CRAN (R 3.6.1)
#  backports              1.1.5     2019-10-02 [1] CRAN (R 3.6.1)
#  Biobase              * 2.46.0    2019-10-29 [2] Bioconductor
#  BiocGenerics         * 0.32.0    2019-10-29 [1] Bioconductor
#  BiocParallel         * 1.20.1    2019-12-21 [1] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.6.1)
#  cli                    2.0.0     2019-12-09 [1] CRAN (R 3.6.1)
#  colorout             * 1.2-2     2019-10-31 [1] Github (jalvesaq/colorout@641ed38)
#  colorspace             1.4-1     2019-03-18 [2] CRAN (R 3.6.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.6.1)
#  DelayedArray         * 0.12.0    2019-10-29 [2] Bioconductor
#  digest                 0.6.23    2019-11-23 [1] CRAN (R 3.6.1)
#  dplyr                * 0.8.3     2019-07-04 [1] CRAN (R 3.6.1)
#  fansi                  0.4.0     2018-10-05 [1] CRAN (R 3.6.1)
#  GenomeInfoDb         * 1.22.0    2019-10-29 [1] Bioconductor
#  GenomeInfoDbData       1.2.2     2019-10-28 [2] Bioconductor
#  GenomicRanges        * 1.38.0    2019-10-29 [1] Bioconductor
#  ggplot2                3.2.1     2019-08-10 [1] CRAN (R 3.6.1)
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
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.6.1)
#  Matrix                 1.2-17    2019-03-22 [3] CRAN (R 3.6.1)
#  matrixStats          * 0.55.0    2019-09-07 [1] CRAN (R 3.6.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.6.1)
#  pillar                 1.4.3     2019-12-20 [1] CRAN (R 3.6.1)
#  pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 3.6.1)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.6.1)
#  promises               1.1.0     2019-10-04 [1] CRAN (R 3.6.1)
#  purrr                  0.3.3     2019-10-18 [2] CRAN (R 3.6.1)
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.6.1)
#  Rcpp                   1.0.3     2019-11-08 [1] CRAN (R 3.6.1)
#  RCurl                  1.95-4.12 2019-03-04 [2] CRAN (R 3.6.1)
#  rlang                  0.4.2     2019-11-23 [1] CRAN (R 3.6.1)
#  rmote                * 0.3.4     2019-10-31 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot              1.3-2     2018-01-03 [2] CRAN (R 3.6.1)
#  S4Vectors            * 0.24.1    2019-12-01 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.6.1)
#  servr                  0.15      2019-08-07 [1] CRAN (R 3.6.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.6.1)
#  SingleCellExperiment * 1.8.0     2019-10-29 [2] Bioconductor
#  SummarizedExperiment * 1.16.1    2019-12-19 [1] Bioconductor
#  tibble                 2.1.3     2019-06-06 [1] CRAN (R 3.6.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.6.1)
#  utf8                   1.1.4     2018-05-24 [1] CRAN (R 3.6.1)
#  vctrs                  0.2.1     2019-12-17 [1] CRAN (R 3.6.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.6.1)
#  xfun                   0.11      2019-11-12 [1] CRAN (R 3.6.1)
#  XVector                0.26.0    2019-10-29 [1] Bioconductor
#  zeallot                0.1.0     2018-01-28 [1] CRAN (R 3.6.1)
#  zlibbioc               1.32.0    2019-10-29 [2] Bioconductor
#
# [1] /users/lcollado/R/3.6.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/library

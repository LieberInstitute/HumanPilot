library('SingleCellExperiment')
library('here')
library('sessioninfo')

## Load data
load(here(
    'Analysis',
    'Human_DLPFC_Visium_processedData_sce_scran.Rdata'
))

## For building the checking function
path <- here('Analysis', 'Layer_Guesses', 'First_Round')
merged_name <- 'Merged'

path <- here('Analysis', 'Layer_Guesses', 'Second_Round')
merged_name <- 'Combined2'


## Read all the csv files
files <- dir(path, pattern = 'csv$', full.names = TRUE)
names(files) <- gsub('.*_|\\.csv', '', dir(path, pattern = 'csv$'))
tabs <- lapply(files, function(f) {
    tab <-
        read.csv(
            f,
            header = TRUE,
            stringsAsFactors = FALSE,
            na.strings = ''
        )
    tab$key <- paste0(tab$sample_name,
        '_',
        tab$spot_name)
    return(tab)
})



table(table(unlist(lapply(tabs[-which(names(tabs) == merged_name)], '[[', 'key'))))
## First round
# 1     2
# 18007  4961

## Second round
# 1
# 24361

## First round continued
table(table(unlist(lapply(tabs[-which(names(tabs) %in% c('151507', merged_name))], '[[', 'key'))))
# 1
# 22968

sapply(tabs, nrow)
# 151507 151669 151673 151508 151670 151674 Merged
# 4961   7857   3611   4381   3484   3635  22968

## For the first round, the data for sample 151507 is duplicated,
## but it's the same
m <- match(tabs[[1]]$key, tabs[[2]]$key)
stopifnot(identical(tabs[[1]]$layer, tabs[[2]]$layer[m]))


## Combine individual slides and check vs merged
merged <- tabs[[which(names(tabs) == merged_name)]]
combined <- do.call(rbind, tabs[-which(names(tabs) == merged_name)])
combined <- combined[!duplicated(combined), ]
rownames(combined) <- NULL

identical(combined, merged) # FALSE =(

m <- match(combined$key, merged$key)
stopifnot(!any(is.na(m)))

identical(combined$key, merged$key[m])


## As a function
check_guesses <- function(path, merged_name) {
    ## Locate and read in the csv files
    files <- dir(path, pattern = 'csv$', full.names = TRUE)
    names(files) <-
        gsub('.*_|\\.csv', '', dir(path, pattern = 'csv$'))
    tabs <- lapply(files, function(f) {
        tab <-
            read.csv(
                f,
                header = TRUE,
                stringsAsFactors = FALSE,
                na.strings = ''
            )
        tab$key <- paste0(tab$sample_name,
            '_',
            tab$spot_name)
        return(tab)
    })
    
    ## Combine individual slides and check vs merged
    merged <- tabs[[which(names(tabs) == merged_name)]]
    combined <-
        do.call(rbind, tabs[-which(names(tabs) == merged_name)])
    combined <- combined[!duplicated(combined), ]
    rownames(combined) <- NULL
    
    ## Re-order
    m <- match(combined$key, merged$key)
    stopifnot(!any(is.na(m)))
    
    ## Check
    identical(combined$key, merged$key[m])
}


check_guesses(
    path = here('Analysis', 'Layer_Guesses', 'First_Round'),
    merged_name = 'Merged'
)
# TRUE

check_guesses(
    path = here('Analysis', 'Layer_Guesses', 'Second_Round'),
    merged_name = 'Combined2'
)
# TRUE

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
# colorspace             1.4-1     2019-03-18 [1] CRAN (R 3.6.1)
# cowplot                1.0.0     2019-07-11 [1] CRAN (R 3.6.1)
# crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.6.1)
# DelayedArray         * 0.12.1    2019-12-17 [1] Bioconductor
# dplyr                  0.8.3     2019-07-04 [1] CRAN (R 3.6.1)
# fansi                  0.4.0     2018-10-05 [1] CRAN (R 3.6.1)
# GenomeInfoDb         * 1.22.0    2019-10-29 [1] Bioconductor
# GenomeInfoDbData       1.2.2     2019-12-18 [1] Bioconductor
# GenomicRanges        * 1.38.0    2019-10-29 [1] Bioconductor
# ggplot2                3.2.1     2019-08-10 [1] CRAN (R 3.6.1)
# glue                   1.3.1     2019-03-12 [1] CRAN (R 3.6.1)
# gtable                 0.3.0     2019-03-25 [1] CRAN (R 3.6.1)
# here                 * 0.1       2017-05-28 [1] CRAN (R 3.6.1)
# IRanges              * 2.20.1    2019-11-20 [1] Bioconductor
# lattice                0.20-38   2018-11-04 [2] CRAN (R 3.6.2)
# lazyeval               0.2.2     2019-03-15 [1] CRAN (R 3.6.1)
# lifecycle              0.1.0     2019-08-01 [1] CRAN (R 3.6.1)
# magrittr               1.5       2014-11-22 [1] CRAN (R 3.6.1)
# Matrix                 1.2-18    2019-11-27 [2] CRAN (R 3.6.2)
# matrixStats          * 0.55.0    2019-09-07 [1] CRAN (R 3.6.1)
# munsell                0.5.0     2018-06-12 [1] CRAN (R 3.6.1)
# packrat                0.5.0     2018-11-14 [1] CRAN (R 3.6.1)
# pillar                 1.4.3     2019-12-20 [1] CRAN (R 3.6.2)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 3.6.1)
# purrr                  0.3.3     2019-10-18 [1] CRAN (R 3.6.1)
# R6                     2.4.1     2019-11-12 [1] CRAN (R 3.6.1)
# Rcpp                   1.0.3     2019-11-08 [1] CRAN (R 3.6.1)
# RCurl                  1.95-4.12 2019-03-04 [1] CRAN (R 3.6.0)
# rlang                  0.4.2     2019-11-23 [1] CRAN (R 3.6.1)
# rprojroot              1.3-2     2018-01-03 [1] CRAN (R 3.6.1)
# rsconnect              0.8.16    2019-12-13 [1] CRAN (R 3.6.1)
# rstudioapi             0.10      2019-03-19 [1] CRAN (R 3.6.1)
# S4Vectors            * 0.24.1    2019-12-01 [1] Bioconductor
# scales                 1.1.0     2019-11-18 [1] CRAN (R 3.6.1)
# sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.6.1)
# SingleCellExperiment * 1.8.0     2019-10-29 [1] Bioconductor
# SummarizedExperiment * 1.16.1    2019-12-20 [1] Bioconductor
# tibble                 2.1.3     2019-06-06 [1] CRAN (R 3.6.1)
# tidyselect             0.2.5     2018-10-11 [1] CRAN (R 3.6.1)
# withr                  2.1.2     2018-03-15 [1] CRAN (R 3.6.1)
# XVector                0.26.0    2019-10-29 [1] Bioconductor
# zlibbioc               1.32.0    2019-10-29 [1] Bioconductor
#
# [1] D:/Documents/R/win-library/3.6
# [2] D:/R/R-3.6.2/library

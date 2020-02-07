library('SingleCellExperiment')
library('here')
library('jaffelab')
library('sessioninfo')

dir.create('pdf', showWarnings = FALSE)
dir.create('rda', showWarnings = FALSE)

## Load data
load(here(
    'Analysis',
    'Human_DLPFC_Visium_processedData_sce_scran.Rdata'
))

## For plotting
source(here('Analysis', 'spatialLIBD_global_plot_code.R'))
genes <- paste0(rowData(sce)$gene_name, '; ', rowData(sce)$gene_id)

genes_to_plot <- c(
    'ATP1A2',
    'FABP7',
    'ADCYAP1',
    'NEFM',
    ## Originally it was NEFN, but that doesn't exist, fixed the typo
    'NEFH',
    'PVALB',
    'VAT1L',
    'NTNG2',
    'CUX2',
    'ENC1',
    'RORB',
    'RXFP1',
    'RPRM',
    'ETV1',
    'PCP4',
    'B3GALT2',
    'CCK',
    'MBP' ## From the figure, not on the Slack message
)

## Load sig_genes data
load('rda/layer_sig_genes.Rdata', verbose = TRUE)

## Looking for NEFN
sort(unique(sig_genes$gene[grep('ne', tolower(sig_genes$gene))]))
# [1] "NECAB1"   "NECAB2"   "NEFH"     "NEFL"     "NEFM"     "NEUROD6"  "SERPINE2"

## After editing genes_to_plot, we are good to go
stopifnot(all(genes_to_plot %in% sig_genes$gene))

## Prepare the data needed for making the plots
sig_genes_sub <- sig_genes[match(genes_to_plot, sig_genes$gene), ]
sig_genes_unique <- splitit(sig_genes_sub$ensembl)

## For the titles
sig_genes_df <- sig_genes_sub
sig_genes_df$in_rows <-
    sapply(sig_genes_df$in_rows, paste0, collapse = ';')
sig_genes_df$results <-
    sapply(sig_genes_df$results, paste0, collapse = ';')


## Make gene grid plots
pdf_dir <- 'pdf/gene_grid/RNAscope'
dir.create(pdf_dir, showWarnings = FALSE, recursive = TRUE)

## Only make the plots for the unique ones
## and only for the last 4 samples
samples_to_plot <- tail(unique(sce$sample_name), 4)
assayname <- 'logcounts'

## From https://gist.githubusercontent.com/mages/5339689/raw/2aaa482dfbbecbfcb726525a3d81661f9d802a8e/add.alpha.R
add.alpha <- function(col, alpha = 1) {
    if (missing(col))
        stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb) / 255, 2,
        function(x)
            rgb(x[1], x[2], x[3], alpha = alpha))
}

for (j in samples_to_plot) {
    # j <-  samples_to_plot[1]
    dir.create(file.path(pdf_dir, j), showWarnings = FALSE)
    
    max_UMI <-
        max(assays(sce)[[assayname]][names(sig_genes_unique), sce$sample_name %in% j])
    x <-
        assays(sce)[[assayname]][names(sig_genes_unique), sce$sample_name %in% j]
    min_UMI <- min(as.vector(x)[as.vector(x) > 0])
    
    for (i in match(names(sig_genes_unique), sig_genes_sub$ensembl)) {
        # i <- 1
        # i <- 15 ## PCP4
        # i <- 11 ## RORB
        message(
            paste(
                Sys.time(),
                'making the plot for',
                i,
                'gene',
                sig_genes_sub$gene[i],
                'minUMI:',
                min_UMI,
                'maxUMI:',
                max_UMI
            )
        )
        
        p <- sce_image_grid_gene(
            sce[, sce$sample_name == j],
            geneid = paste0(sig_genes_sub$gene[i], '; ', sig_genes_sub$ensembl[i]),
            return_plots = TRUE,
            ... = gsub('top', 'r', gsub(
                'Layer', 'L', sig_genes_df$results[i]
            )),
            spatial = TRUE,
            assayname = assayname,
            minCount = 0
        )
        
        p2 <- p[[1]] + scale_fill_gradientn(
            colors = c('aquamarine4', 'springgreen', 'goldenrod', 'red'),
            na.value = add.alpha('black', 0.175),
            name = assayname,
            ## Scales code borrowed from Comparison_SpatialDE_genes_pooled.html
            values = scales::rescale(c(min_UMI, 2, 4, max_UMI))
        ) +
            scale_color_gradientn(
                colors =  c('aquamarine4', 'springgreen', 'goldenrod', 'red'),
                na.value =  add.alpha('black', 0.175),
                name = assayname,
                values = scales::rescale(c(min_UMI, 2, 4, max_UMI))
            ) + theme(legend.title = element_text(size = 20),
                legend.text = element_text(size = 15))
        pdf(
            file.path(pdf_dir,
                j,
                paste0(
                    sig_genes_sub$gene[i],
                    gsub('top', 'r', gsub(
                        'Layer', 'L', sig_genes_df$results[i]
                    )),
                    '.pdf'
                )),
            useDingbats = FALSE,
            height = 8,
            width = 9.5
        )
        print(p2)
        dev.off()
    }
}


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
#  date     2020-02-07
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
#  googledrive            1.0.0     2019-08-19 [1] CRAN (R 3.6.1)
#  gtable                 0.3.0     2019-03-25 [2] CRAN (R 3.6.1)
#  here                 * 0.1       2017-05-28 [1] CRAN (R 3.6.1)
#  htmltools              0.4.0     2019-10-04 [1] CRAN (R 3.6.1)
#  htmlwidgets            1.5.1     2019-10-08 [1] CRAN (R 3.6.1)
#  httpuv                 1.5.2     2019-09-11 [1] CRAN (R 3.6.1)
#  IRanges              * 2.20.1    2019-11-20 [1] Bioconductor
#  jaffelab             * 0.99.29   2019-11-04 [1] Github (LieberInstitute/jaffelab@a7d87cb)
#  jsonlite               1.6       2018-12-07 [2] CRAN (R 3.6.1)
#  labeling               0.3       2014-08-23 [2] CRAN (R 3.6.1)
#  later                  1.0.0     2019-10-04 [1] CRAN (R 3.6.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.6.1)
#  lazyeval               0.2.2     2019-03-15 [2] CRAN (R 3.6.1)
#  limma                  3.42.0    2019-10-29 [1] Bioconductor
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
#  rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 3.6.1)
#  RColorBrewer           1.1-2     2014-12-07 [2] CRAN (R 3.6.1)
#  Rcpp                   1.0.3     2019-11-08 [1] CRAN (R 3.6.1)
#  RCurl                  1.95-4.12 2019-03-04 [2] CRAN (R 3.6.1)
#  rlang                  0.4.2     2019-11-23 [1] CRAN (R 3.6.1)
#  rmote                * 0.3.4     2019-10-31 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot              1.3-2     2018-01-03 [2] CRAN (R 3.6.1)
#  S4Vectors            * 0.24.1    2019-12-01 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.6.1)
#  scatterplot3d          0.3-41    2018-03-14 [1] CRAN (R 3.6.1)
#  segmented              1.0-0     2019-06-17 [2] CRAN (R 3.6.1)
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

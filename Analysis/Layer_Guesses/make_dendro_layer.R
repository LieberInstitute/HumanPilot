library('SingleCellExperiment')
library('here')
library('readxl')
library('Polychrome')
library('rafalib')
library('sessioninfo')

## Functions derived from this script, to make it easier to resume the work
sce_layer_file <-
    here('Analysis', 'Layer_Guesses', 'rda', 'sce_layer.Rdata')
if (file.exists(sce_layer_file))
    load(sce_layer_file, verbose = TRUE)
source(here('Analysis', 'Layer_Guesses', 'layer_specificity_functions.R'),
    echo = TRUE)


d <- dist(t(assays(sce_layer)$logcounts))
h <- hclust(d)

labs <- paste0(
    sce_layer$subject,
    '.',
    sce_layer$position,
    '.',
    ifelse(sce_layer$replicate == 1, 'A', 'B')
)

pdf('pdf/dendro_layer.pdf',
    width = 22,
    useDingbats = FALSE)
pal <- RColorBrewer::brewer.pal(12, 'Paired')
pal[11:12] <- c('maroon1', 'magenta3')
palette(pal)
par(mar = c(1.5, 5, 2, 1) + 0.1)
myplclust(
    h,
    labels = gsub('ayer', '', sce_layer$layer_guess),
    lab.col = as.numeric(as.factor(labs)),
    main = 'Layers by sample (brain ID, µm position, replicate)',
    cex = 2,
    cex.axis = 2,
    cex.lab = 2,
    cex.main = 2
)
legend(
    'topright',
    unique(labs),
    col = pal,
    lwd = 4,
    bty = 'n',
    cex = 1.5,
    ncol = 3
)
dev.off()

pdf('pdf/dendro_layer_extra.pdf',
    width = 22,
    useDingbats = FALSE)
palette(Polychrome::palette36.colors(7))
par(mar = c(8.5, 5, 2, 1) + 0.1)
myplclust(
    h,
    labels = labs,
    lab.col = as.numeric(sce_layer$layer_guess),
    main = 'Layers by sample (brain ID, µm position, replicate)',
    cex = 2,
    cex.axis = 2,
    cex.lab = 2,
    cex.main = 2
)
legend(
    'topright',
    c(paste0('L', 1:6), 'WM'),
    col = Polychrome::palette36.colors(7)[c(2:7, 1)],
    lwd = 4,
    bty = 'n',
    cex = 2,
    ncol = 3
)
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
#  date     2020-02-24
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
#  magrittr               1.5      2014-11-22 [1] CRAN (R 3.6.1)
#  Matrix                 1.2-17   2019-03-22 [3] CRAN (R 3.6.1)
#  matrixStats          * 0.55.0   2019-09-07 [1] CRAN (R 3.6.1)
#  munsell                0.5.0    2018-06-12 [2] CRAN (R 3.6.1)
#  pillar                 1.4.3    2019-12-20 [1] CRAN (R 3.6.1)
#  pkgconfig              2.0.3    2019-09-22 [1] CRAN (R 3.6.1)
#  png                    0.1-7    2013-12-03 [2] CRAN (R 3.6.1)
#  Polychrome           * 1.2.3    2019-08-01 [1] CRAN (R 3.6.1)
#  promises               1.1.0    2019-10-04 [1] CRAN (R 3.6.1)
#  R6                     2.4.1    2019-11-12 [2] CRAN (R 3.6.1)
#  rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 3.6.1)
#  RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 3.6.1)
#  Rcpp                   1.0.3    2019-11-08 [1] CRAN (R 3.6.1)
#  RCurl                  1.98-1.1 2020-01-19 [2] CRAN (R 3.6.1)
#  readxl               * 1.3.1    2019-03-13 [2] CRAN (R 3.6.1)
#  rlang                  0.4.2    2019-11-23 [1] CRAN (R 3.6.1)
#  rmote                * 0.3.4    2019-10-31 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot              1.3-2    2018-01-03 [2] CRAN (R 3.6.1)
#  S4Vectors            * 0.24.1   2019-12-01 [1] Bioconductor
#  scales                 1.1.0    2019-11-18 [2] CRAN (R 3.6.1)
#  scatterplot3d          0.3-41   2018-03-14 [1] CRAN (R 3.6.1)
#  servr                  0.15     2019-08-07 [1] CRAN (R 3.6.1)
#  sessioninfo          * 1.1.1    2018-11-05 [1] CRAN (R 3.6.1)
#  SingleCellExperiment * 1.8.0    2019-10-29 [2] Bioconductor
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

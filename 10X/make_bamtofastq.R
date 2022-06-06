library('sgejobs')
library('sessioninfo')

dirs <- dir(pattern = '^1')
stopifnot(length(dirs) == 12)

job_loop(
    loops = list(sample = dirs),
    name = 'bamtofastq',
    cores = 4,
    queue = 'bluejay',
    memory = '20G',
    create_shell = TRUE,
    logdir = 'logs_bamtofastq'
)
dir.create('logs_bamtofastq', showWarnings = FALSE)

# /dcl02/lieber/ajaffe/SpatialTranscriptomics/bamtofastq-1.2.0 --nthreads=4 /dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/10X/${sample}/${sample}_mRNA.bam /dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/10X/${sample}/fastq


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 3.6.1 Patched (2019-09-06 r77160)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2019-12-03
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package     * version    date       lib source
#  assertthat    0.2.1      2019-03-21 [2] CRAN (R 3.6.1)
#  backports     1.1.5      2019-10-02 [1] CRAN (R 3.6.1)
#  cli           1.1.0      2019-03-19 [1] CRAN (R 3.6.1)
#  codetools     0.2-16     2018-12-24 [3] CRAN (R 3.6.1)
#  colorout    * 1.2-2      2019-09-26 [1] Github (jalvesaq/colorout@641ed38)
#  colorspace    1.4-1      2019-03-18 [2] CRAN (R 3.6.1)
#  crayon        1.3.4      2017-09-16 [1] CRAN (R 3.6.1)
#  curl          4.2        2019-09-24 [1] CRAN (R 3.6.1)
#  digest        0.6.22     2019-10-21 [1] CRAN (R 3.6.1)
#  dplyr         0.8.3      2019-07-04 [1] CRAN (R 3.6.1)
#  ggplot2       3.2.1      2019-08-10 [1] CRAN (R 3.6.1)
#  glue          1.3.1      2019-03-12 [1] CRAN (R 3.6.1)
#  gtable        0.3.0      2019-03-25 [2] CRAN (R 3.6.1)
#  hms           0.5.2      2019-10-30 [1] CRAN (R 3.6.1)
#  htmltools     0.4.0      2019-10-04 [1] CRAN (R 3.6.1)
#  htmlwidgets   1.5.1      2019-10-08 [1] CRAN (R 3.6.1)
#  httpuv        1.5.2      2019-09-11 [1] CRAN (R 3.6.1)
#  jsonlite      1.6        2018-12-07 [2] CRAN (R 3.6.1)
#  later         1.0.0      2019-10-04 [1] CRAN (R 3.6.1)
#  lattice       0.20-38    2018-11-04 [3] CRAN (R 3.6.1)
#  lazyeval      0.2.2      2019-03-15 [2] CRAN (R 3.6.1)
#  lifecycle     0.1.0      2019-08-01 [1] CRAN (R 3.6.1)
#  lubridate     1.7.4      2018-04-11 [1] CRAN (R 3.6.1)
#  magrittr      1.5        2014-11-22 [1] CRAN (R 3.6.1)
#  munsell       0.5.0      2018-06-12 [2] CRAN (R 3.6.1)
#  pillar        1.4.2      2019-06-29 [1] CRAN (R 3.6.1)
#  pkgconfig     2.0.3      2019-09-22 [1] CRAN (R 3.6.1)
#  png           0.1-7      2013-12-03 [2] CRAN (R 3.6.1)
#  promises      1.1.0      2019-10-04 [1] CRAN (R 3.6.1)
#  pryr          0.1.4      2018-02-18 [2] CRAN (R 3.6.1)
#  purrr         0.3.3      2019-10-18 [1] CRAN (R 3.6.1)
#  R6            2.4.0      2019-02-14 [2] CRAN (R 3.6.1)
#  Rcpp          1.0.3      2019-11-08 [1] CRAN (R 3.6.1)
#  readr         1.3.1      2018-12-21 [1] CRAN (R 3.6.1)
#  remotes       2.1.0.9000 2019-09-26 [1] Github (MangoTheCat/remotes@2f1a040)
#  rlang         0.4.1      2019-10-24 [1] CRAN (R 3.6.1)
#  rmote       * 0.3.4      2019-09-26 [1] Github (cloudyr/rmote@fbce611)
#  scales        1.0.0      2018-08-09 [2] CRAN (R 3.6.1)
#  servr         0.15       2019-08-07 [1] CRAN (R 3.6.1)
#  sessioninfo * 1.1.1      2018-11-05 [1] CRAN (R 3.6.1)
#  sgejobs     * 0.99.1     2019-11-07 [1] Github (LieberInstitute/sgejobs@f5ab0ca)
#  stringi       1.4.3      2019-03-12 [2] CRAN (R 3.6.1)
#  stringr       1.4.0      2019-02-10 [1] CRAN (R 3.6.1)
#  tibble        2.1.3      2019-06-06 [1] CRAN (R 3.6.1)
#  tidyr         1.0.0      2019-09-11 [1] CRAN (R 3.6.1)
#  tidyselect    0.2.5      2018-10-11 [2] CRAN (R 3.6.1)
#  vctrs         0.2.0      2019-07-05 [1] CRAN (R 3.6.1)
#  withr         2.1.2      2018-03-15 [2] CRAN (R 3.6.1)
#  xfun          0.10       2019-10-01 [1] CRAN (R 3.6.1)
#  zeallot       0.1.0      2018-01-28 [1] CRAN (R 3.6.1)
#
# [1] /users/lcollado/R/3.6
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6/R/3.6/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6/R/3.6/lib64/R/library
# >

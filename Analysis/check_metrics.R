###
library(jaffelab)
metricFiles = list.files(
    "/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/10X",
    pattern = "metrics_summary_csv.csv",
    full = TRUE,
    recur = TRUE
)
names(metricFiles) = ss(metricFiles, "/", 8)

metrics = sapply(metricFiles, read.csv, as.is = TRUE)


## with high mean rates of exonic alignments (mean: XX%, IQR: XX-XX%)
summary(as.numeric(gsub('\\%', '',
    unlist(metrics['Reads.Mapped.Confidently.to.Exonic.Regions', ]))))
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 81.50   82.50   83.55   83.31   84.28   85.00


## Add info about cell segmentation to the metrics table
## Code originally from
## https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/Layer_Guesses/spots_per_layer.R
library('SingleCellExperiment')
library('here')
library('dplyr')
library('sessioninfo')

## Load data
load(here(
    'Analysis',
    'Human_DLPFC_Visium_processedData_sce_scran.Rdata'
))

summary(sce$cell_count)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.000   1.000   3.000   3.281   4.000  27.000


mean_cells <-
    group_by(as.data.frame(colData(sce)), sample_name) %>% summarize(mean_cells = mean(cell_count))
## Check the order
stopifnot(identical(colnames(metrics), as.character(mean_cells$sample_name)))
metrics <-
    rbind(metrics, 'Mean.Cells.Per.Spot' = mean_cells$mean_cells)

## Add info about cells per layer
load(here('Analysis', 'Layer_Guesses', 'rda', 'layer_guess_tab.Rdata'),
    verbose = TRUE)

layer_guess_tab$cell_count = sce$cell_count[match(layer_guess_tab$key, sce$key)]

prop_sample <- group_by(layer_guess_tab, sample_name) %>%
    summarize(prop0 = mean(cell_count == 0),
        prop1 = mean(cell_count == 1))

metrics <-
    rbind(metrics, 'Proportion.0.Cells.Per.Spot' = prop_sample$prop0)
metrics <-
    rbind(metrics, 'Proportion.1.Cell.Per.Spot' = prop_sample$prop1)

## Add some more sample info
stopifnot(identical(colnames(metrics), as.character(unique(sce$sample_name))))
m <- match(unique(sce$sample_name), sce$sample_name)

metrics <- rbind(metrics, 'Brain.Number' = sce$subject[m])
metrics <- rbind(metrics, 'Position' = paste(sce$position[m], 'µm'))
metrics <- rbind(metrics, 'Replicate' = sce$replicate[m])

library('LIBDpheno')
m_pheno <- match(sce$subject[m], toxicant[['2019-09-03']]$brnum)
metrics <-
    rbind(metrics, 'Age.Death' = toxicant[['2019-09-03']]$agedeath[m_pheno])
metrics <-
    rbind(metrics, 'Sex' = as.character(toxicant[['2019-09-03']]$sex[m_pheno]))
metrics <-
    rbind(metrics, 'Primary.Diagnosis' = as.character(toxicant[['2019-09-03']]$primarydx[m_pheno]))


metrics <- cbind('Variable' = rownames(metrics), metrics)
write.table(metrics, sep = "\t",
    file = "visium_dlpfc_pilot_sample_metrics.tsv", row.names = FALSE)

## We sequenced each sample to a mean depth of XX reads
summary(as.numeric(gsub(',', '', unlist(metrics['Number.of.Reads', -1])))) / 1e6
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 229.3   269.3   291.1   346.1   327.7   822.2

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
#  date     2020-02-17                                 
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source                                    
#  assertthat             0.2.1     2019-03-21 [2] CRAN (R 3.6.1)                            
#  backports              1.1.5     2019-10-02 [1] CRAN (R 3.6.1)                            
#  Biobase              * 2.46.0    2019-10-29 [2] Bioconductor                              
#  BiocGenerics         * 0.32.0    2019-10-29 [1] Bioconductor                              
#  BiocParallel         * 1.20.1    2019-12-21 [1] Bioconductor                              
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.6.1)                            
#  callr                  3.4.0     2019-12-09 [1] CRAN (R 3.6.1)                            
#  cellranger             1.1.0     2016-07-27 [1] CRAN (R 3.6.1)                            
#  cli                    2.0.0     2019-12-09 [1] CRAN (R 3.6.1)                            
#  colorout             * 1.2-2     2019-10-31 [1] Github (jalvesaq/colorout@641ed38)        
#  colorspace             1.4-1     2019-03-18 [2] CRAN (R 3.6.1)                            
#  cowplot              * 1.0.0     2019-07-11 [1] CRAN (R 3.6.1)                            
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.6.1)                            
#  curl                   4.3       2019-12-02 [1] CRAN (R 3.6.1)                            
#  data.table             1.12.8    2019-12-09 [1] CRAN (R 3.6.1)                            
#  DelayedArray         * 0.12.0    2019-10-29 [2] Bioconductor                              
#  desc                   1.2.0     2018-05-01 [2] CRAN (R 3.6.1)                            
#  devtools               2.2.1     2019-09-24 [1] CRAN (R 3.6.1)                            
#  digest                 0.6.23    2019-11-23 [1] CRAN (R 3.6.1)                            
#  dplyr                * 0.8.3     2019-07-04 [1] CRAN (R 3.6.1)                            
#  DT                     0.11      2019-12-19 [1] CRAN (R 3.6.1)                            
#  ellipsis               0.3.0     2019-09-20 [1] CRAN (R 3.6.1)                            
#  fansi                  0.4.0     2018-10-05 [1] CRAN (R 3.6.1)                            
#  fastmap                1.0.1     2019-10-08 [1] CRAN (R 3.6.1)                            
#  forcats                0.4.0     2019-02-17 [1] CRAN (R 3.6.1)                            
#  foreign                0.8-72    2019-08-02 [3] CRAN (R 3.6.1)                            
#  fs                     1.3.1     2019-05-06 [1] CRAN (R 3.6.1)                            
#  gargle                 0.4.0     2019-10-04 [1] CRAN (R 3.6.1)                            
#  GenomeInfoDb         * 1.22.0    2019-10-29 [1] Bioconductor                              
#  GenomeInfoDbData       1.2.2     2019-10-28 [2] Bioconductor                              
#  GenomicRanges        * 1.38.0    2019-10-29 [1] Bioconductor                              
#  ggplot2              * 3.2.1     2019-08-10 [1] CRAN (R 3.6.1)                            
#  glue                   1.3.1     2019-03-12 [1] CRAN (R 3.6.1)                            
#  googleAuthR            1.1.1     2019-09-09 [1] CRAN (R 3.6.1)                            
#  googledrive            1.0.0     2019-08-19 [1] CRAN (R 3.6.1)                            
#  gtable                 0.3.0     2019-03-25 [2] CRAN (R 3.6.1)                            
#  haven                  2.2.0     2019-11-08 [1] CRAN (R 3.6.1)                            
#  here                 * 0.1       2017-05-28 [1] CRAN (R 3.6.1)                            
#  hms                    0.5.2     2019-10-30 [2] CRAN (R 3.6.1)                            
#  htmltools              0.4.0     2019-10-04 [1] CRAN (R 3.6.1)                            
#  htmlwidgets            1.5.1     2019-10-08 [1] CRAN (R 3.6.1)                            
#  httpuv                 1.5.2     2019-09-11 [1] CRAN (R 3.6.1)                            
#  httr                   1.4.1     2019-08-05 [1] CRAN (R 3.6.1)                            
#  IRanges              * 2.20.1    2019-11-20 [1] Bioconductor                              
#  jaffelab             * 0.99.29   2019-11-04 [1] Github (LieberInstitute/jaffelab@a7d87cb) 
#  jsonlite               1.6       2018-12-07 [2] CRAN (R 3.6.1)                            
#  later                  1.0.0     2019-10-04 [1] CRAN (R 3.6.1)                            
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.6.1)                            
#  lazyeval               0.2.2     2019-03-15 [2] CRAN (R 3.6.1)                            
#  LIBDpheno            * 0.99.69   2019-11-04 [1] Github (LieberInstitute/LIBDpheno@01dda5e)
#  limma                  3.42.0    2019-10-29 [1] Bioconductor                              
#  lmtest                 0.9-37    2019-04-30 [2] CRAN (R 3.6.1)                            
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.6.1)                            
#  MASS                   7.3-51.4  2019-03-31 [3] CRAN (R 3.6.1)                            
#  Matrix                 1.2-17    2019-03-22 [3] CRAN (R 3.6.1)                            
#  matrixStats          * 0.55.0    2019-09-07 [1] CRAN (R 3.6.1)                            
#  memoise                1.1.0     2017-04-21 [2] CRAN (R 3.6.1)                            
#  mime                   0.8       2019-12-19 [1] CRAN (R 3.6.1)                            
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.6.1)                            
#  openxlsx               4.1.4     2019-12-06 [1] CRAN (R 3.6.1)                            
#  pillar                 1.4.3     2019-12-20 [1] CRAN (R 3.6.1)                            
#  pkgbuild               1.0.6     2019-10-09 [2] CRAN (R 3.6.1)                            
#  pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 3.6.1)                            
#  pkgload                1.0.2     2018-10-29 [2] CRAN (R 3.6.1)                            
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.6.1)                            
#  Polychrome           * 1.2.3     2019-08-01 [1] CRAN (R 3.6.1)                            
#  prettyunits            1.0.2     2015-07-13 [1] CRAN (R 3.6.1)                            
#  processx               3.4.1     2019-07-18 [1] CRAN (R 3.6.1)                            
#  promises               1.1.0     2019-10-04 [1] CRAN (R 3.6.1)                            
#  ps                     1.3.0     2018-12-21 [2] CRAN (R 3.6.1)                            
#  purrr                  0.3.3     2019-10-18 [2] CRAN (R 3.6.1)                            
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.6.1)                            
#  rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 3.6.1)                            
#  RColorBrewer           1.1-2     2014-12-07 [2] CRAN (R 3.6.1)                            
#  Rcpp                   1.0.3     2019-11-08 [1] CRAN (R 3.6.1)                            
#  RCurl                  1.95-4.12 2019-03-04 [2] CRAN (R 3.6.1)                            
#  readxl                 1.3.1     2019-03-13 [2] CRAN (R 3.6.1)                            
#  remotes                2.1.0     2019-06-24 [1] CRAN (R 3.6.1)                            
#  rio                    0.5.16    2018-11-26 [1] CRAN (R 3.6.1)                            
#  rlang                  0.4.2     2019-11-23 [1] CRAN (R 3.6.1)                            
#  rmote                * 0.3.4     2019-10-31 [1] Github (cloudyr/rmote@fbce611)            
#  rprojroot              1.3-2     2018-01-03 [2] CRAN (R 3.6.1)                            
#  S4Vectors            * 0.24.1    2019-12-01 [1] Bioconductor                              
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.6.1)                            
#  scatterplot3d          0.3-41    2018-03-14 [1] CRAN (R 3.6.1)                            
#  segmented              1.0-0     2019-06-17 [2] CRAN (R 3.6.1)                            
#  servr                  0.15      2019-08-07 [1] CRAN (R 3.6.1)                            
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.6.1)                            
#  shiny                  1.4.0     2019-10-10 [1] CRAN (R 3.6.1)                            
#  shinycsv               0.99.8    2019-11-04 [1] Github (LieberInstitute/shinycsv@7f5e49d) 
#  SingleCellExperiment * 1.8.0     2019-10-29 [2] Bioconductor                              
#  stringi                1.4.3     2019-03-12 [2] CRAN (R 3.6.1)                            
#  SummarizedExperiment * 1.16.1    2019-12-19 [1] Bioconductor                              
#  testthat               2.3.1     2019-12-01 [1] CRAN (R 3.6.1)                            
#  tibble                 2.1.3     2019-06-06 [1] CRAN (R 3.6.1)                            
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.6.1)                            
#  usethis                1.5.1     2019-07-04 [1] CRAN (R 3.6.1)                            
#  utf8                   1.1.4     2018-05-24 [1] CRAN (R 3.6.1)                            
#  vcd                    1.4-4     2017-12-06 [1] CRAN (R 3.6.1)                            
#  vctrs                  0.2.1     2019-12-17 [1] CRAN (R 3.6.1)                            
#  viridisLite          * 0.3.0     2018-02-01 [2] CRAN (R 3.6.1)                            
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.6.1)                            
#  xfun                   0.11      2019-11-12 [1] CRAN (R 3.6.1)                            
#  xtable                 1.8-4     2019-04-21 [2] CRAN (R 3.6.1)                            
#  XVector                0.26.0    2019-10-29 [1] Bioconductor                              
#  zeallot                0.1.0     2018-01-28 [1] CRAN (R 3.6.1)                            
#  zip                    2.0.4     2019-09-01 [1] CRAN (R 3.6.1)                            
#  zlibbioc               1.32.0    2019-10-29 [2] Bioconductor                              
#  zoo                    1.8-6     2019-05-28 [2] CRAN (R 3.6.1)                            
# 
# [1] /users/lcollado/R/3.6.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/library

library('SingleCellExperiment')
library('broom')
library('ggplot2')
library('cowplot')
library('sessioninfo')

## Load the sce object from sce_scran.R
load('Human_DLPFC_Visium_processedData_sce_scran.Rdata',
    verbose = TRUE)
dir.create('pdf', showWarnings = FALSE, mode = '0770')

##
d <- as.data.frame(colData(sce))



fit <- lm(cell_count ~ sum_umi + sample_name, data = d)
summary(fit)
# Call:
# lm(formula = cell_count ~ sum_umi + sample_name, data = d)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -4.3426 -1.8301 -0.5018  1.1205 23.8984
#
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)
# (Intercept) -3.252e+02  2.469e+01  -13.17   <2e-16 ***
# sum_umi      1.200e-04  6.279e-06   19.11   <2e-16 ***
# sample_name  2.164e-03  1.629e-04   13.29   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 2.603 on 47678 degrees of freedom
# Multiple R-squared:  0.018,     Adjusted R-squared:  0.01795
# F-statistic: 436.8 on 2 and 47678 DF,  p-value: < 2.2e-16

fit2 <- lm(cell_count ~ sum_gene + sample_name, data = d)
summary(fit2)
# Call:
# lm(formula = cell_count ~ sum_gene + sample_name, data = d)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -4.4418 -1.8187 -0.5343  1.1184 24.0868
#
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)
# (Intercept) -2.661e+02  2.452e+01  -10.85   <2e-16 ***
# sum_gene     4.273e-04  1.647e-05   25.94   <2e-16 ***
# sample_name  1.772e-03  1.618e-04   10.95   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 2.595 on 47678 degrees of freedom
# Multiple R-squared:  0.02424,   Adjusted R-squared:  0.0242
# F-statistic: 592.2 on 2 and 47678 DF,  p-value: < 2.2e-16

fit3 <- lm(cell_count ~ sum_umi + sum_gene + sample_name, data = d)
summary(fit3)

# Call:
# lm(formula = cell_count ~ sum_umi + sum_gene + sample_name, data = d)
#
# Residuals:
#     Min      1Q  Median      3Q     Max
# -4.6083 -1.7697 -0.5134  1.1255 24.3671
#
# Coefficients:
#               Estimate Std. Error t value Pr(>|t|)
# (Intercept) -3.362e+02  2.431e+01  -13.83   <2e-16 ***
# sum_umi     -1.130e-03  3.309e-05  -34.16   <2e-16 ***
# sum_gene     3.350e-03  8.709e-05   38.46   <2e-16 ***
# sample_name  2.226e-03  1.604e-04   13.88   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
# Residual standard error: 2.564 on 47677 degrees of freedom
# Multiple R-squared:  0.04755,   Adjusted R-squared:  0.04749
# F-statistic: 793.4 on 3 and 47677 DF,  p-value: < 2.2e-16


anova(fit, fit2, fit3)
# Analysis of Variance Table
#
# Model 1: cell_count ~ sum_umi + sum_gene + sample_name
# Model 2: cell_count ~ sum_gene + sample_name
# Model 3: cell_count ~ sum_umi + sum_gene + sample_name
#   Res.Df    RSS Df Sum of Sq      F    Pr(>F)
# 1  47677 313343
# 2  47678 321012 -1   -7668.8 1166.9 < 2.2e-16 ***
# 3  47677 313343  1    7668.8 1166.9 < 2.2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

fits <- lapply(split(d, d$sample_name), function(x) {
    lm(cell_count ~ sum_umi + sum_gene, data = x)
})


fits_tab <-
    do.call(rbind, mapply(
        function(x, y)
            cbind(tidy(x), sample_name = y),
        fits,
        names(fits),
        SIMPLIFY = FALSE
    ))
rownames(fits_tab) <- NULL

data.frame(
    'sample_name' = names(fits),
    'sum_umi' = p.adjust(subset(fits_tab, term == 'sum_umi')$p.value, 'fdr') < 0.05,
    'sum_gene' = p.adjust(subset(fits_tab, term == 'sum_gene')$p.value, 'fdr') < 0.05
)
#    sample_name sum_umi sum_gene
# 1       151507    TRUE     TRUE
# 2       151508    TRUE     TRUE
# 3       151509    TRUE     TRUE
# 4       151510    TRUE     TRUE
# 5       151669    TRUE     TRUE
# 6       151670    TRUE     TRUE
# 7       151671   FALSE     TRUE
# 8       151672   FALSE    FALSE
# 9       151673    TRUE     TRUE
# 10      151674    TRUE     TRUE
# 11      151675    TRUE     TRUE
# 12      151676    TRUE     TRUE





## Should switch the outcome to the more continuous variable

# fit <- lm(log(sum_umi) ~ cell_count + sample_name, data = d)
# summary(fit)

fits <- lapply(split(d, d$sample_name), function(x) {
    models <- list(
        'umi_vs_cell' = lm(log(sum_umi) ~ cell_count, data = x),
        'gene_vs_cell' = lm(log(sum_gene) ~ cell_count, data = x),
        'gene_vs_umi' = lm(log(sum_gene) ~ log(sum_umi), data = x)
    )
})


fits_tab <-
    do.call(rbind, mapply(function(x, y) {
        res <- do.call(rbind, lapply(x, function(z)
            tidy(z)[2, ]))
        res$sample_name <- y
        res$model <- names(x)
        return(res)
    },
        fits,
        names(fits),
        SIMPLIFY = FALSE))
rownames(fits_tab) <- NULL

fits_tab_l <- split(fits_tab, fits_tab$model)
fits_tab_l <- lapply(fits_tab_l, function(x) {
    x$FDR_raw <- p.adjust(x$p.value, 'fdr')
    x$FDR <- format(x$FDR_raw, format = 'e', digits = 2)
    x$signif <-
        factor(x$FDR_raw < 0.05, levels = c('TRUE', 'FALSE'))
    return(x)
})

p <-
    ggplot(d, aes(x = cell_count, y = log(sum_umi))) + geom_point() +
    geom_smooth(method = 'lm') + facet_grid( ~ sample_name) +
    theme_bw(base_size = 20) +
    geom_text(
        data = fits_tab_l[['umi_vs_cell']],
        mapping = aes(
            x = -Inf,
            y = -Inf,
            label = FDR,
            color = signif
        ),
        hjust = -0.2,
        vjust = -1,
        size = 5
    )

p2 <-
    ggplot(d, aes(x = cell_count, y = log(sum_gene))) + geom_point() +
    geom_smooth(method = 'lm') + facet_grid(~ sample_name) +
    theme_bw(base_size = 20) +
    geom_text(
        data = fits_tab_l[['gene_vs_cell']],
        mapping = aes(
            x = -Inf,
            y = -Inf,
            label = FDR,
            color = signif
        ),
        hjust = -0.2,
        vjust = -1,
        size = 5
    )

p3 <-
    ggplot(d, aes(x = log(sum_umi), y = log(sum_gene))) + geom_point() +
    geom_smooth(method = 'lm') + facet_grid(~ sample_name) +
    theme_bw(base_size = 20) +
    geom_text(
        data = fits_tab_l[['gene_vs_umi']],
        mapping = aes(
            x = -Inf,
            y = -Inf,
            label = FDR,
            color = signif
        ),
        hjust = -0.2,
        vjust = -1,
        size = 5
    )


pdf(
    'pdf/cell_count_vs_umi.pdf',
    width = 40,
    height = 21,
    useDingbats = FALSE
)
plot_grid(p, p2, p3, ncol = 1)
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
#  date     2019-12-11
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  assertthat             0.2.1     2019-03-21 [2] CRAN (R 3.6.1)
#  backports              1.1.5     2019-10-02 [1] CRAN (R 3.6.1)
#  Biobase              * 2.46.0    2019-10-29 [2] Bioconductor
#  BiocGenerics         * 0.32.0    2019-10-29 [1] Bioconductor
#  BiocParallel         * 1.20.0    2019-10-30 [1] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.6.1)
#  broom                * 0.5.2     2019-04-07 [1] CRAN (R 3.6.1)
#  cli                    1.1.0     2019-03-19 [1] CRAN (R 3.6.1)
#  colorout             * 1.2-2     2019-10-31 [1] Github (jalvesaq/colorout@641ed38)
#  colorspace             1.4-1     2019-03-18 [2] CRAN (R 3.6.1)
#  cowplot              * 1.0.0     2019-07-11 [1] CRAN (R 3.6.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.6.1)
#  DelayedArray         * 0.12.0    2019-10-29 [2] Bioconductor
#  digest                 0.6.22    2019-10-21 [1] CRAN (R 3.6.1)
#  dplyr                  0.8.3     2019-07-04 [1] CRAN (R 3.6.1)
#  generics               0.0.2     2018-11-29 [1] CRAN (R 3.6.1)
#  GenomeInfoDb         * 1.22.0    2019-10-29 [1] Bioconductor
#  GenomeInfoDbData       1.2.2     2019-10-28 [2] Bioconductor
#  GenomicRanges        * 1.38.0    2019-10-29 [1] Bioconductor
#  ggplot2              * 3.2.1     2019-08-10 [1] CRAN (R 3.6.1)
#  glue                   1.3.1     2019-03-12 [1] CRAN (R 3.6.1)
#  gtable                 0.3.0     2019-03-25 [2] CRAN (R 3.6.1)
#  htmltools              0.4.0     2019-10-04 [1] CRAN (R 3.6.1)
#  htmlwidgets            1.5.1     2019-10-08 [1] CRAN (R 3.6.1)
#  httpuv                 1.5.2     2019-09-11 [1] CRAN (R 3.6.1)
#  IRanges              * 2.20.0    2019-10-29 [1] Bioconductor
#  jsonlite               1.6       2018-12-07 [2] CRAN (R 3.6.1)
#  later                  1.0.0     2019-10-04 [1] CRAN (R 3.6.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.6.1)
#  lazyeval               0.2.2     2019-03-15 [2] CRAN (R 3.6.1)
#  lifecycle              0.1.0     2019-08-01 [1] CRAN (R 3.6.1)
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.6.1)
#  Matrix                 1.2-17    2019-03-22 [3] CRAN (R 3.6.1)
#  matrixStats          * 0.55.0    2019-09-07 [1] CRAN (R 3.6.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.6.1)
#  nlme                   3.1-141   2019-08-01 [3] CRAN (R 3.6.1)
#  pillar                 1.4.2     2019-06-29 [1] CRAN (R 3.6.1)
#  pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 3.6.1)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.6.1)
#  promises               1.1.0     2019-10-04 [1] CRAN (R 3.6.1)
#  purrr                  0.3.3     2019-10-18 [2] CRAN (R 3.6.1)
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.6.1)
#  Rcpp                   1.0.3     2019-11-08 [1] CRAN (R 3.6.1)
#  RCurl                  1.95-4.12 2019-03-04 [2] CRAN (R 3.6.1)
#  rlang                  0.4.1     2019-10-24 [1] CRAN (R 3.6.1)
#  rmote                * 0.3.4     2019-10-31 [1] Github (cloudyr/rmote@fbce611)
#  S4Vectors            * 0.24.0    2019-10-29 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.6.1)
#  servr                  0.15      2019-08-07 [1] CRAN (R 3.6.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.6.1)
#  SingleCellExperiment * 1.8.0     2019-10-29 [2] Bioconductor
#  SummarizedExperiment * 1.16.0    2019-10-29 [1] Bioconductor
#  tibble                 2.1.3     2019-06-06 [1] CRAN (R 3.6.1)
#  tidyr                  1.0.0     2019-09-11 [1] CRAN (R 3.6.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.6.1)
#  vctrs                  0.2.0     2019-07-05 [1] CRAN (R 3.6.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.6.1)
#  xfun                   0.10      2019-10-01 [1] CRAN (R 3.6.1)
#  XVector                0.26.0    2019-10-29 [1] Bioconductor
#  zeallot                0.1.0     2018-01-28 [1] CRAN (R 3.6.1)
#  zlibbioc               1.32.0    2019-10-29 [2] Bioconductor
#
# [1] /users/lcollado/R/3.6.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/library

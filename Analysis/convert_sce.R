# module load conda_R/3.6.x
library('SingleCellExperiment')
library('ggplot2')
library('sessioninfo')

## load rse list
load("Human_DLPFC_Visium_processedData_rseList.rda", verbose = TRUE)

sceList <- lapply(rseList, function(rse) {
    SingleCellExperiment(
        assays = list(counts = assays(rse)$umis),
        rowData = rowData(rse),
        colData = colData(rse),
        metadata = metadata(rse)
    )
})
save(sceList, file = 'Human_DLPFC_Visium_processedData_sceList.Rdata')

sce <- do.call(cbind, sceList)
metadata(sce) <- list('image' = do.call(rbind, metadata(sce)))
table(colData(sce)$sample_name)
# 151507 151508 151509 151510 151669 151670 151671 151672 151673 151674 151675
#   4226   4384   4789   4634   3661   3498   4110   4015   3639   3673   3592
# 151676
#   3460

## Add design info
study <-
    read.csv(
        '/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/image_index_10xID.csv'
    )

## same order
stopifnot(identical(match(names(sceList), study$X10xID),
    seq_len(length(sceList))))


gsub('JHU_|_rep.*', '', study$Description)
gsub('.*position | µm', '', study$Description)
sce$subject <- rep(gsub('JHU_|_rep.*', '', study$Description),
    table(colData(sce)$sample_name))
with(colData(sce), table(sample_name, subject))
#            subject
# sample_name Br5292 Br5595 Br8100
#      151507   4226      0      0
#      151508   4384      0      0
#      151509   4789      0      0
#      151510   4634      0      0
#      151669      0   3661      0
#      151670      0   3498      0
#      151671      0   4110      0
#      151672      0   4015      0
#      151673      0      0   3639
#      151674      0      0   3673
#      151675      0      0   3592
#      151676      0      0   3460

sce$position <- rep(gsub('.*position | µm', '', study$Description),
    table(colData(sce)$sample_name))
with(colData(sce), table(sample_name, position))
#            position
# sample_name    0  300
#      151507 4226    0
#      151508 4384    0
#      151509    0 4789
#      151510    0 4634
#      151669 3661    0
#      151670 3498    0
#      151671    0 4110
#      151672    0 4015
#      151673 3639    0
#      151674 3673    0
#      151675    0 3592
#      151676    0 3460

sce$replicate <- rep(gsub('.*0', '', study$Rep),
    table(colData(sce)$sample_name))
with(colData(sce), table(sample_name, replicate))
#            replicate
# sample_name    1    2
#      151507 4226    0
#      151508    0 4384
#      151509 4789    0
#      151510    0 4634
#      151669 3661    0
#      151670    0 3498
#      151671 4110    0
#      151672    0 4015
#      151673 3639    0
#      151674    0 3673
#      151675 3592    0
#      151676    0 3460

## For blocking later with scran
sce$subject_position <- paste0(sce$subject, '_pos', sce$position)

save(sce, file = 'Human_DLPFC_Visium_processedData_sce.Rdata')


## Plotting function
#' @param ... Parameters passed to paste0() for the plot title
sce_image_clus <- function(sce,
    sampleid,
    clustervar,
    colors = c(
        "#b2df8a",
        "#e41a1c",
        "#377eb8",
        "#4daf4a",
        "#ff7f00",
        "gold",
        "#a65628",
        "#999999",
        "black",
        "grey",
        "white",
        "purple"
    ),
    ...) {
    d <- as.data.frame(colData(sce[, sce$sample_name == sampleid]))
    p <- ggplot(d,
        aes(
            x = imagecol,
            y = imagerow,
            fill = factor(!!sym(clustervar))
        )) +
        geom_spatial(
            data = subset(metadata(sce)$image, sample == sampleid),
            aes(grob = grob),
            x = 0.5,
            y = 0.5
        ) +
        geom_point(shape = 21,
            size = 1.25,
            stroke = 0.25) +
        coord_cartesian(expand = FALSE) +
        scale_fill_manual(values = colors) +
        xlim(0, max(sce$width)) +
        ylim(max(sce$height), 0) +
        xlab("") + ylab("") +
        labs(fill = clustervar) +
        guides(fill = guide_legend(override.aes = list(size = 3))) +
        ggtitle(paste0(sampleid, ...)) +
        theme_set(theme_bw(base_size = 10)) +
        theme(
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text = element_blank(),
            axis.ticks = element_blank()
        )
    return(p)
}

pdf('test_plot.pdf')
sce_image_clus(sce, '151673', 'Cluster', ... = ' k 6')
dev.off()

unlink('test_plot.pdf')


sort_clusters <- function(clusters, map_subset = NULL) {
    if (is.null(map_subset)) {
        map_subset <- rep(TRUE, length(clusters))
    }
    map <-
        rank(length(clusters[map_subset]) - table(clusters[map_subset]), ties.method = 'first')
    res <- map[clusters]
    factor(res)
}

# > clusters <- c(1, 1, 1, 2, 2, 3, 4, 4, 4, 4)
# > sort_clusters(clusters)
# 1 1 1 2 2 3 4 4 4 4
# 2 2 2 3 3 4 1 1 1 1
# Levels: 1 2 3 4
# > sort_clusters(sort_clusters(clusters))
# 2 2 2 3 3 4 1 1 1 1
# 2 2 2 3 3 4 1 1 1 1
# Levels: 1 2 3 4


sce_image_grid <-
    function(sce,
        clusters,
        pdf_file,
        sort_clust = TRUE,
        colors = NULL,
        return_plots = FALSE,
        ...) {
        if (is.null(colors)) {
            ## Original ones
            # colors <- c("#b2df8a","#e41a1c","#377eb8","#4daf4a","#ff7f00","gold",
            # "#a65628", "#999999", "black", "grey", "white", "purple")
            
            ## From https://medialab.github.io/iwanthue/
            ## which I found the link to from
            ## https://stackoverflow.com/questions/15282580/how-to-generate-a-number-of-most-distinctive-colors-in-r
            
            ## Used the colorblind friendly and default palette
            n_clus <- length(unique(clusters))
            ## Color-blind friendly
            # colors <- if(n_clus > 10) c("#573487", "#96b43d", "#5771dd", "#d39830", "#a874d7", "#64c36f", "#c86ebd", "#47bb8a", "#892862", "#33d4d1", "#db5358", "#6a98e4", "#c55e32", "#516bba", "#b3ad52", "#d55f90", "#588234", "#b8465e", "#a97a35", "#a44739") else if (n_clus <= 6) c("#98a441", "#6778d0", "#50b47b","#9750a1", "#ba6437", "#ba496b") else c("#bfa239", "#5a3a8e", "#799f44", "#c771c4", "#54b06c", "#b0457b", "#43c9b0", "#b8434e", "#6e81da", "#b86738")
            ## Default one
            # colors <- if(n_clus > 10) c("#aab539", "#6f66d7", "#59b94d", "#bd54c1", "#598124", "#d1478f", "#58bd91", "#d4465a", "#48bbd2", "#d6542d", "#6288cc", "#d69938", "#8d61ab", "#a3b26a", "#da8dc6", "#407e4a", "#a24b66", "#86712e", "#e39178", "#a75932") else if (n_clus <= 6) c("#b88f40", "#7a75cd", "#6ca74d", "#c45ca2", "#49adaa", "#cb584c") else c("#6ab64c", "#8761cc", "#c1a942", "#688bcc", "#d35238", "#4db598", "#c361aa", "#677e39", "#c9566e", "#c07b44")
            
            ## Default one if n_clus > 12, otherwise original colors (re-ordered a bit)
            colors <-
                if (n_clus > 12)
                    c(
                        "#aab539",
                        "#6f66d7",
                        "#59b94d",
                        "#bd54c1",
                        "#598124",
                        "#d1478f",
                        "#58bd91",
                        "#d4465a",
                        "#48bbd2",
                        "#d6542d",
                        "#6288cc",
                        "#d69938",
                        "#8d61ab",
                        "#a3b26a",
                        "#da8dc6",
                        "#407e4a",
                        "#a24b66",
                        "#86712e",
                        "#e39178",
                        "#a75932"
                    )
            else
                c(
                    "#377eb8",
                    "gold",
                    "#ff7f00",
                    "#e41a1c",
                    "#4daf4a",
                    "#b2df8a",
                    "#a65628",
                    "#999999",
                    "black",
                    "grey",
                    "white",
                    "purple"
                )
            names(colors) <- seq_len(length(colors))
            
        }
        sce$Clus <- if (sort_clust)
            sort_clusters(clusters)
        else
            clusters
        plots <- lapply(unique(sce$sample_name), function(sampleid) {
            sce_image_clus(sce, sampleid, 'Clus', colors = colors, ...)
        })
        if (!return_plots) {
            pdf(pdf_file, height = 24, width = 36)
            print(cowplot::plot_grid(plotlist = plots))
            dev.off()
            return(pdf_file)
        }
        else {
            return(plots)
        }
    }

sce_image_grid_comp <-
    function(sce,
        clus1,
        clus2,
        pdf_file,
        map_subset = NULL,
        colors = NULL,
        return_plots = FALSE,
        ...) {
        clus1 <- sort_clusters(clus1, map_subset)
        clus2 <- sort_clusters(clus2, map_subset)
        if (is.null(colors)) {
            colors <- c('FALSE' = 'red', 'TRUE' = 'grey80')
        }
        clus_agree <-
            factor(clus1 == clus2, levels = c('FALSE', 'TRUE'))
        sce_image_grid(
            sce,
            clus_agree,
            pdf_file,
            sort_clust = FALSE,
            colors = colors,
            return_plots = return_plots,
            ...
        )
    }


sce_image_grid_by_clus <-
    function(sce, clusters, pdf_file, colors = NULL, ...) {
        if (is.null(colors)) {
            colors <- c('FALSE' = 'transparent', 'TRUE' = 'red')
        }
        clusters_uni <- sort(unique(clusters))
        pdf(pdf_file, height = 24, width = 36)
        lapply(clusters_uni, function(clus) {
            curr_clus <- factor(clusters == clus, levels = c('FALSE', 'TRUE'))
            plots <-
                sce_image_grid(
                    sce,
                    curr_clus,
                    sort_clust = FALSE,
                    colors = colors,
                    return_plots = TRUE,
                    ... = paste(..., '- cluster', clus)
                )
            print(cowplot::plot_grid(plotlist = plots))
            return(clus)
        })
        dev.off()
        return(pdf_file)
    }

## Plotting layer and function
save(geom_spatial,
    sce_image_clus,
    sort_clusters,
    sce_image_grid,
    file = 'geom_spatial.Rdata')

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
#  date     2019-11-11
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  assertthat             0.2.1     2019-03-21 [2] CRAN (R 3.6.1)
#  backports              1.1.5     2019-10-02 [1] CRAN (R 3.6.1)
#  Biobase              * 2.46.0    2019-10-29 [2] Bioconductor
#  BiocGenerics         * 0.32.0    2019-10-29 [1] Bioconductor
#  BiocParallel         * 1.20.0    2019-10-30 [1] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.6.1)
#  cli                    1.1.0     2019-03-19 [1] CRAN (R 3.6.1)
#  colorout             * 1.2-2     2019-10-31 [1] Github (jalvesaq/colorout@641ed38)
#  colorspace             1.4-1     2019-03-18 [2] CRAN (R 3.6.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.6.1)
#  DelayedArray         * 0.12.0    2019-10-29 [2] Bioconductor
#  digest                 0.6.22    2019-10-21 [1] CRAN (R 3.6.1)
#  dplyr                  0.8.3     2019-07-04 [1] CRAN (R 3.6.1)
#  fansi                  0.4.0     2018-10-05 [1] CRAN (R 3.6.1)
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
#  labeling               0.3       2014-08-23 [2] CRAN (R 3.6.1)
#  later                  1.0.0     2019-10-04 [1] CRAN (R 3.6.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.6.1)
#  lazyeval               0.2.2     2019-03-15 [2] CRAN (R 3.6.1)
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.6.1)
#  Matrix                 1.2-17    2019-03-22 [3] CRAN (R 3.6.1)
#  matrixStats          * 0.55.0    2019-09-07 [1] CRAN (R 3.6.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.6.1)
#  pillar                 1.4.2     2019-06-29 [1] CRAN (R 3.6.1)
#  pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 3.6.1)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.6.1)
#  promises               1.1.0     2019-10-04 [1] CRAN (R 3.6.1)
#  purrr                  0.3.3     2019-10-18 [2] CRAN (R 3.6.1)
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.6.1)
#  Rcpp                   1.0.2     2019-07-25 [1] CRAN (R 3.6.1)
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
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.6.1)
#  utf8                   1.1.4     2018-05-24 [1] CRAN (R 3.6.1)
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

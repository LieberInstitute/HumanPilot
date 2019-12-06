# screen -S scran
# qrsh -l mem_free=60G,h_vmem=60G,h_fsize=100G -pe local 4
# module load conda_R/3.6.x
library('SingleCellExperiment')
library('scran')
library('scater')
library('BiocParallel')
library('PCAtools')
library('igraph')
library('ggplot2')
library('cowplot')
library('jaffelab') ## for ss(), splitit(), myplclust(); actually they are from rafalib
library('pheatmap')
library('sessioninfo')

# library('zinbwave')
# library('clusterExperiment')
#
# library('RColorBrewer')

dir.create('pdf_scran', showWarnings = FALSE)
dir.create('rda_scran', showWarnings = FALSE)

## From convert_sce.R
load('geom_spatial.Rdata', verbose = TRUE)

## For resuming and other analyses
if (!file.exists('Human_DLPFC_Visium_processedData_sce_scran.Rdata')) {
    load('Human_DLPFC_Visium_processedData_sce.Rdata', verbose = TRUE)
    ## From
    ## http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#2_setting_up_the_data
    qcstats <- perCellQCMetrics(sce)
    qcfilter <- quickPerCellQC(qcstats)
    colSums(as.matrix(qcfilter))
    # low_lib_size low_n_features        discard
    #          451            510            534
    
    with(qcfilter, table(low_lib_size, low_n_features, discard))
    # , , discard = FALSE
    #
    #             low_n_features
    # low_lib_size FALSE  TRUE
    #        FALSE 47147     0
    #        TRUE      0     0
    #
    # , , discard = TRUE
    #
    #             low_n_features
    # low_lib_size FALSE  TRUE
    #        FALSE     0    83
    #        TRUE     24   427
    
    table(sce$sample_name[qcfilter$discard])
    # 151507 151508 151509 151510 151669 151670 151671 151672 151673 151674 151675
    #     45    111     59     38     28     30     59    139      7      3      9
    # 151676
    #      6
    
    ## Plot discarded umis
    load('geom_spatial.Rdata', verbose = TRUE)
    sce$discard <- qcfilter$discard
    plots_discard <-
        lapply(unique(sce$sample_name), function(sampleid) {
            sce_image_clus(sce, sampleid, 'discard', colors = c('light blue', 'red'))
        })
    pdf('pdf_scran/discarded_cells_grid.pdf',
        height = 24,
        width = 36)
    plot_grid(plotlist = plots_discard)
    dev.off()
    
    ## We have decided not to filter umis since they seem to be layer specific
    # sce <- sce[,!qcfilter$discard]
    # summary(qcfilter$discard)
    #    Mode   FALSE    TRUE
    # logical   47147     534
    
    ## From
    ## http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#3_normalizing_cell-specific_biases
    set.seed(20191112)
    Sys.time()
    clusters <- quickCluster(
        sce,
        BPPARAM = MulticoreParam(4),
        block = sce$subject_position,
        block.BPPARAM = MulticoreParam(4)
    )
    Sys.time()
    ## Takes about 2 minutes
    # [1] "2019-11-13 10:56:34 EST"
    # [1] "2019-11-13 10:57:56 EST"
    
    
    sce <-
        computeSumFactors(sce, clusters = clusters, BPPARAM = MulticoreParam(4))
    Sys.time()
    ## Takes about 3 minutes
    # [1] "2019-11-13 11:00:32 EST"
    
    summary(sizeFactors(sce))
    #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # 0.1329  0.5613  0.8552  1.0000  1.2695  6.7378
    sce <- logNormCounts(sce)
    
    ## From
    ## http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#4_variance_modelling
    dec <- modelGeneVar(sce,
        block = sce$subject_position,
        BPPARAM = MulticoreParam(4))
    Sys.time()
    ## Takes about 30 secs to get here from the computeSumFactors() step
    # [1] "2019-11-13 11:00:59 EST"
    
    pdf('pdf_scran/modelGeneVar.pdf', useDingbats = FALSE)
    mapply(function(block, blockname) {
        plot(
            block$mean,
            block$total,
            xlab = "Mean log-expression",
            ylab = "Variance",
            main = blockname
        )
        curve(metadata(block)$trend(x),
            col = "blue",
            add = TRUE)
    }, dec$per.block, names(dec$per.block))
    dev.off()
    
    top.hvgs <- getTopHVGs(dec, prop = 0.1)
    length(top.hvgs)
    # [1] 1942
    
    ## Basically the same
    # top.hvgs2 <- getTopHVGs(dec, n=2000)
    # table(top.hvgs %in% top.hvgs2)
    
    top.hvgs.fdr5 <- getTopHVGs(dec, fdr.threshold = 0.05)
    length(top.hvgs.fdr5)
    # [1] 13842
    
    top.hvgs.fdr1 <- getTopHVGs(dec, fdr.threshold = 0.01)
    length(top.hvgs.fdr1)
    # [1] 12393
    
    ## FDR-based selection returns too many genes
    
    set.seed(20191112)
    Sys.time()
    sce <- runPCA(sce, subset_row = top.hvgs)
    Sys.time()
    ## Takes about 2 minutes
    # [1] "2019-11-13 11:02:24 EST"
    # [1] "2019-11-13 11:04:33 EST"
    
    reducedDimNames(sce)
    
    ## PCs don't have sd = 1
    summary(apply(reducedDim(sce, 'PCA'), 2, sd))
    #   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    # 0.8810  0.8948  0.9143  1.1404  1.0950  4.2064
    
    summary(apply(reducedDim(sce, 'PCA'), 2, sd))
    #       PC1       PC2       PC3       PC4       PC5       PC6       PC7       PC8
    # 4.2063825 2.9181994 2.4577468 1.9652432 1.7231614 1.4718300 1.3936265 1.3330686
    #       PC9      PC10      PC11      PC12      PC13      PC14      PC15      PC16
    # 1.1943737 1.1733142 1.1550694 1.1171480 1.1067906 1.0595489 1.0156686 0.9805711
    #      PC17      PC18      PC19      PC20      PC21      PC22      PC23      PC24
    # 0.9676315 0.9423317 0.9346976 0.9278445 0.9266340 0.9222307 0.9197716 0.9182687
    #      PC25      PC26      PC27      PC28      PC29      PC30      PC31      PC32
    # 0.9147040 0.9139619 0.9123189 0.9092247 0.9078313 0.9029707 0.9022757 0.9009750
    #      PC33      PC34      PC35      PC36      PC37      PC38      PC39      PC40
    # 0.8999425 0.8980520 0.8973509 0.8957044 0.8952632 0.8947088 0.8933552 0.8930551
    #      PC41      PC42      PC43      PC44      PC45      PC46      PC47      PC48
    # 0.8915526 0.8897350 0.8890296 0.8872520 0.8871032 0.8855887 0.8832843 0.8824588
    #      PC49      PC50
    # 0.8821567 0.8810461
    
    ## Means are 0 though
    summary(colMeans(reducedDim(sce, 'PCA')))
    #       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
    # -1.225e-14 -3.569e-15  1.417e-15  2.647e-15  8.394e-15  2.893e-14
    
    ## 2019-12-06 edits: add tsne, umap and spot cell numbers
    
    
    ## From https://github.com/davismcc/scater/blob/master/R/runTSNE.R#L85
    ## I see that the default perplexity will be 50
    # > mat <- scater:::.get_mat_from_sce(sce, exprs_values = 'logcounts', dimred = 'PCA', n_dimred = NULL)
    # > dim(mat)
    # [1] 47681    50
    # > min(50, floor(nrow(mat) / 5))
    # [1] 50
    Sys.time()
    set.seed(20191206)
    sce <- runTSNE(sce, dimred = 'PCA', name = 'TSNE_perplexity50', perplexity = 50)
    Sys.time()
    ## Takes about 14 min
    # [1] "2019-12-06 14:07:53 EST"
    # [1] "2019-12-06 14:21:55 EST"
    
    Sys.time()
    set.seed(20191206)
    sce <- runTSNE(sce, dimred = 'PCA', name = 'TSNE_perplexity5', perplexity = 5)
    Sys.time()
    ## Takes about 10 min
    # [1] "2019-12-06 14:22:15 EST"
    # [1] "2019-12-06 14:32:35 EST"
    
    Sys.time()
    set.seed(20191206)
    sce <- runTSNE(sce, dimred = 'PCA', name = 'TSNE_perplexity20', perplexity = 20)
    Sys.time()
    ## Takes about 12 min
    # [1] "2019-12-06 14:44:38 EST"
    
    Sys.time()
    set.seed(20191206)
    sce <- runTSNE(sce, dimred = 'PCA', name = 'TSNE_perplexity80', perplexity = 80)
    Sys.time()
    
    ## Takes about 15 min
    # [1] "2019-12-06 14:44:59 EST"
    # [1] "2019-12-06 15:00:45 EST"
    
    ## From https://github.com/davismcc/scater/blob/master/R/runUMAP.R#L65
    ## looks like the default n_neighbors is 15
    Sys.time()
    set.seed(20191206)
    sce <- runUMAP(sce, dimred = 'PCA', name = 'UMAP_neighbors15')
    Sys.time()
    
    ## Takes about 2 mins
    # [1] "2019-12-06 15:07:14 EST"
    # [1] "2019-12-06 15:08:59 EST"
    
    
    ## Read in the number of cells per spot
    cells <- do.call(rbind, lapply(dir('Histology'), function(sampleid) {
        x <- read.csv(file.path('Histology', sampleid, 'tissue_spot_counts.csv'))
        x$key <- paste0(sampleid, '_', x$barcode)
        return(x[, c('key', 'count')])
    }))
    
    ## Used in plotly code in spatialLIBD
    sce$key <- paste0(sce$sample_name, '_', colnames(sce))
    m <- match(sce$key, cells$key)
    stopifnot(!all(is.na(m)))
    sce$cell_count <- cells$count[m]
    
    tapply(sce$cell_count, sce$sample_name, summary)
    # $`151507`
    #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    #   0.000   1.000   2.000   2.202   3.000  19.000
    #
    # $`151508`
    #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    #    0.00    1.00    3.00    3.31    5.00   24.00
    #
    # $`151509`
    #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    #   0.000   1.000   3.000   2.994   4.000  21.000
    #
    # $`151510`
    #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    #   0.000   2.000   3.000   3.202   5.000  13.000
    #
    # $`151669`
    #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    #   0.000   2.000   2.000   2.633   4.000  10.000
    #
    # $`151670`
    #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    #   0.000   3.000   5.000   5.595   8.000  21.000
    #
    # $`151671`
    #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    #   0.000   2.000   3.000   2.733   4.000  13.000
    #
    # $`151672`
    #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    #   0.000   1.000   2.000   2.066   3.000  13.000
    #
    # $`151673`
    #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    #   0.000   3.000   4.000   4.531   6.000  27.000
    #
    # $`151674`
    #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    #   0.000   2.000   3.000   3.971   5.000  22.000
    #
    # $`151675`
    #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    #   0.000   2.000   3.000   3.457   5.000  24.000
    #
    # $`151676`
    #    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
    #   0.000   1.000   2.000   3.243   4.000  20.000
    save(sce, top.hvgs, file = 'Human_DLPFC_Visium_processedData_sce_scran.Rdata')
    
} else {
    load('Human_DLPFC_Visium_processedData_sce_scran.Rdata',
        verbose = TRUE)
}

## From
## https://osca.bioconductor.org/dimensionality-reduction.html#using-the-elbow-point
percent.var <- attr(reducedDim(sce), "percentVar")
chosen.elbow <- PCAtools::findElbowPoint(percent.var)
chosen.elbow
# [1] 2

pdf('pdf_scran/PCA_var_explained.pdf', useDingbats = FALSE)
plot(percent.var, xlab = "PC", ylab = "Variance explained (%)")
abline(v = chosen.elbow, col = "red")
dev.off()



## From
## http://bioconductor.org/packages/release/bioc/vignettes/scran/inst/doc/scran.html#5_automated_pc_choice
set.seed(20191112)
Sys.time()
sced <- denoisePCA(sce, dec, subset.row = top.hvgs)
Sys.time()
ncol(reducedDim(sced, "PCA"))
# [1] 9

## Takes about a minute to run
# [1] "2019-11-13 11:16:25 EST"
# [1] "2019-11-13 11:17:20 EST"



# Sys.time()
# choices <- getClusteredPCs(reducedDim(sce))
# Sys.time()

## This takes forever to run... killed after about 20 hours
# [1] "2019-11-12 12:53:43 EST"
## Killed around:
# [1] "2019-11-13 10:56:34 EST"

# npcs <- metadata(choices)$chosen
# npcs
# reducedDim(sce, "PCAsub") <- reducedDim(sce, "PCA")[, seq_len(npcs), drop = FALSE]
#
#
# pdf('pdf_scran/PC_choices.pdf', useDingbats = FALSE)
# plot(choices$n.pcs, choices$n.clusters,
#     xlab="Number of PCs", ylab="Number of clusters")
# abline(a=1, b=1, col="red")
# abline(v=metadata(choices)$chosen, col="grey80", lty=2)
# dev.off()

Sys.time()
g_k10 <- buildSNNGraph(sce, k = 10, use.dimred = 'PCA')
Sys.time()

## Takes 2 min
# [1] "2019-11-13 11:21:29 EST"
# [1] "2019-11-13 11:23:22 EST"

Sys.time()
g_walk_k10 <- igraph::cluster_walktrap(g_k10)
clust_k10 <- g_walk_k10$membership
Sys.time()

## This takes longer (about 1 hour)
# [1] "2019-11-13 12:21:06 EST"

clust_k10 <- sort_clusters(g_walk_k10$membership)
save(g_k10, g_walk_k10, file = 'rda_scran/g_k10.Rdata')

table(clust_k10)
# clust_k10
#    1    2    3    4    5    6    7    8    9   10   11   12   13
# 8982 7674 7615 6423 5513 3831 3081 1375 1321 1011  772   50   33

options(width = 200)
addmargins(table(clust_k10, sce$subject_position))
# clust_k10 Br5292_pos0 Br5292_pos300 Br5595_pos0 Br5595_pos300 Br8100_pos0 Br8100_pos300   Sum
#       1          1237          1247         568           711        2904          2315  8982
#       2            69            64        3841          3526          33           141  7674
#       3          3546          3928          15            26          28            72  7615
#       4           245           310         116           125        2776          2851  6423
#       5            45            19        2308          3066          26            49  5513
#       6          1619          2005          91            19          31            66  3831
#       7          1269          1325         169            22         120           176  3081
#       8           385           207          17           601          67            98  1375
#       9            25            62          17            13         545           659  1321
#       10          131            90           0             4         487           299  1011
#       11           19           149           0             0         286           318   772
#       12           11             8          13             7           7             4    50
#       13            9             9           4             5           2             4    33
#       Sum        8610          9423        7159          8125        7312          7052 47681

addmargins(table(clust_k10, sce$subject))
# clust_k10 Br5292 Br5595 Br8100   Sum
#       1     2484   1279   5219  8982
#       2      133   7367    174  7674
#       3     7474     41    100  7615
#       4      555    241   5627  6423
#       5       64   5374     75  5513
#       6     3624    110     97  3831
#       7     2594    191    296  3081
#       8      592    618    165  1375
#       9       87     30   1204  1321
#       10     221      4    786  1011
#       11     168      0    604   772
#       12      19     20     11    50
#       13      18      9      6    33
#       Sum  18033  15284  14364 47681


## Needs 28 colors! (for K = 50 further below, use the same colors here then)
## From https://medialab.github.io/iwanthue/ with the default preset
cols <- c(
    "#de84b0",
    "#78bb40",
    "#9c45bd",
    "#46c06f",
    "#cb3e97",
    "#488733",
    "#b978e9",
    "#beae36",
    "#5c6ade",
    "#db9443",
    "#5985dc",
    "#cf4e32",
    "#43c4c4",
    "#d84068",
    "#5fb88e",
    "#e471d1",
    "#327e58",
    "#7454b1",
    "#a4b266",
    "#964f95",
    "#72722a",
    "#c18cd3",
    "#a06332",
    "#54a4d6",
    "#dc8074",
    "#5465a4",
    "#9f4765",
    "#a09cdf"
)
names(cols) <- seq_len(length(cols))

sce_image_grid(sce, clust_k10, 'pdf_scran/grid_SNN_k10_noXY.pdf', colors = cols)


## Try with another k
Sys.time()
g_k50 <- buildSNNGraph(sce, k = 50, use.dimred = 'PCA')
Sys.time()
## About 12 minutes
# [1] "2019-11-13 15:20:32 EST"
# [1] "2019-11-13 15:31:50 EST"

Sys.time()
g_walk_k50 <- igraph::cluster_walktrap(g_k50)
Sys.time()

## About 1 hour? Nope, closer to a day
# [1] "2019-11-13 15:31:50 EST"
# [1] "2019-11-14 12:05:23 EST"

clust_k50 <- sort_clusters(g_walk_k50$membership)
save(g_k50, g_walk_k50, file = 'rda_scran/g_k50.Rdata')

sce_image_grid(sce, clust_k50, 'pdf_scran/grid_SNN_k50_noXY.pdf', colors = cols)

table(clust_k50, useNA = 'ifany')
# clust_k50
#    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16
# 6020 4750 3926 3764 3726 3704 3243 2284 2278 2259 1677 1676 1341 1228  761  579
#   17   18   19   20   21   22   23   24   25   26   27   28
#  555  497  481  429  426  424  413  386  337  305  142   70

addmargins(table(clust_k50, sce$subject))

# clust_k50 Br5292 Br5595 Br8100   Sum
#       1       20   5964     36  6020
#       2      101   4420    229  4750
#       3     1634    117   2175  3926
#       4      304     56   3404  3764
#       5     3683     13     30  3726
#       6      338   2617    749  3704
#       7      730    333   2180  3243
#       8     2279      0      5  2284
#       9     2033    134    111  2278
#       10     594     40   1625  2259
#       11    1668      6      3  1677
#       12    1217    264    195  1676
#       13     905    157    279  1341
#       14     907    144    177  1228
#       15     408    119    234   761
#       16      89      3    487   579
#       17     251    254     50   555
#       18      56     66    375   497
#       19      68      0    413   481
#       20      42     10    377   429
#       21      70      0    356   426
#       22     100      0    324   424
#       23     135    273      5   413
#       24     142    160     84   386
#       25      69      3    265   337
#       26      41     92    172   305
#       27      90     36     16   142
#       28      59      3      8    70
#       Sum  18033  15284  14364 47681

## For it to be k = 7
clust_k50_k7 <- sort_clusters(igraph::cut_at(g_walk_k50, n = 7))
sce_image_grid(sce, clust_k50_k7, 'pdf_scran/grid_SNN_k50_k7_noXY.pdf')
addmargins(table(clust_k50_k7, sce$subject))
# clust_k50_k7 Br5292 Br5595 Br8100   Sum
#          1    12109    829   9617 22555
#          2      459  13001   1014 14474
#          3     2493    708    962  4163
#          4     2033    134    111  2278
#          5      327      3   1580  1910
#          6      167     79   1017  1263
#          7      445    530     63  1038
#          Sum  18033  15284  14364 47681


addmargins(table(clust_k50, clust_k50_k7))
#          clust_k50_k7
# clust_k50     1     2     3     4     5     6     7   Sum
#       1       0  6020     0     0     0     0     0  6020
#       2       0  4750     0     0     0     0     0  4750
#       3    3926     0     0     0     0     0     0  3926
#       4    3764     0     0     0     0     0     0  3764
#       5    3726     0     0     0     0     0     0  3726
#       6       0  3704     0     0     0     0     0  3704
#       7    3243     0     0     0     0     0     0  3243
#       8    2284     0     0     0     0     0     0  2284
#       9       0     0     0  2278     0     0     0  2278
#       10   2259     0     0     0     0     0     0  2259
#       11   1677     0     0     0     0     0     0  1677
#       12   1676     0     0     0     0     0     0  1676
#       13      0     0  1341     0     0     0     0  1341
#       14      0     0  1228     0     0     0     0  1228
#       15      0     0   761     0     0     0     0   761
#       16      0     0     0     0   579     0     0   579
#       17      0     0     0     0     0     0   555   555
#       18      0     0     0     0     0   497     0   497
#       19      0     0     0     0   481     0     0   481
#       20      0     0     0     0     0   429     0   429
#       21      0     0     0     0   426     0     0   426
#       22      0     0     0     0   424     0     0   424
#       23      0     0     0     0     0     0   413   413
#       24      0     0   386     0     0     0     0   386
#       25      0     0     0     0     0   337     0   337
#       26      0     0   305     0     0     0     0   305
#       27      0     0   142     0     0     0     0   142
#       28      0     0     0     0     0     0    70    70
#       Sum 22555 14474  4163  2278  1910  1263  1038 47681


k50_summ <-
    as.data.frame(
        table(
            'TheirCluster' = sce$Cluster,
            'ClusterK50' = clust_k50,
            'ClusterK50_Cut7' = clust_k50_k7,
            'sample_name' = sce$sample_name
        )
    )
dim(k50_summ)
# [1] 21168     5
k50_summ <- subset(k50_summ, Freq != 0)
dim(k50_summ)
# [1] 905   5
write.csv(k50_summ,
    file = 'rda_scran/k50_summ.csv',
    row.names = FALSE,
    quote = FALSE)

clust_k50_k14 <- sort_clusters(igraph::cut_at(g_walk_k50, n = 14))
sce_image_grid(sce,
    clust_k50_k14,
    'pdf_scran/grid_SNN_k50_k14_noXY.pdf',
    colors = cols)
addmargins(table(clust_k50_k14, sce$subject))
# clust_k50_k14 Br5292 Br5595 Br8100   Sum
#           1      121  10384    265 10770
#           2     5894    603   2383  8880
#           3     5317    130   2205  7652
#           4      898     96   5029  6023
#           5      338   2617    749  3704
#           6     1812    301    456  2569
#           7     2033    134    111  2278
#           8      681    407    506  1594
#           9      157      3    900  1060
#           10     170      0    680   850
#           11     111     13    642   766
#           12     310    257     58   625
#           13      56     66    375   497
#           14     135    273      5   413
#           Sum  18033  15284  14364 47681
# addmargins(table(clust_k50_k14, sce$sample_name))


addmargins(table(clust_k50, clust_k50_k14))
#          clust_k50_k14
# clust_k50     1     2     3     4     5     6     7     8     9    10    11    12    13    14   Sum
#       1    6020     0     0     0     0     0     0     0     0     0     0     0     0     0  6020
#       2    4750     0     0     0     0     0     0     0     0     0     0     0     0     0  4750
#       3       0     0  3926     0     0     0     0     0     0     0     0     0     0     0  3926
#       4       0     0     0  3764     0     0     0     0     0     0     0     0     0     0  3764
#       5       0     0  3726     0     0     0     0     0     0     0     0     0     0     0  3726
#       6       0     0     0     0  3704     0     0     0     0     0     0     0     0     0  3704
#       7       0  3243     0     0     0     0     0     0     0     0     0     0     0     0  3243
#       8       0  2284     0     0     0     0     0     0     0     0     0     0     0     0  2284
#       9       0     0     0     0     0     0  2278     0     0     0     0     0     0     0  2278
#       10      0     0     0  2259     0     0     0     0     0     0     0     0     0     0  2259
#       11      0  1677     0     0     0     0     0     0     0     0     0     0     0     0  1677
#       12      0  1676     0     0     0     0     0     0     0     0     0     0     0     0  1676
#       13      0     0     0     0     0  1341     0     0     0     0     0     0     0     0  1341
#       14      0     0     0     0     0  1228     0     0     0     0     0     0     0     0  1228
#       15      0     0     0     0     0     0     0   761     0     0     0     0     0     0   761
#       16      0     0     0     0     0     0     0     0   579     0     0     0     0     0   579
#       17      0     0     0     0     0     0     0     0     0     0     0   555     0     0   555
#       18      0     0     0     0     0     0     0     0     0     0     0     0   497     0   497
#       19      0     0     0     0     0     0     0     0   481     0     0     0     0     0   481
#       20      0     0     0     0     0     0     0     0     0     0   429     0     0     0   429
#       21      0     0     0     0     0     0     0     0     0   426     0     0     0     0   426
#       22      0     0     0     0     0     0     0     0     0   424     0     0     0     0   424
#       23      0     0     0     0     0     0     0     0     0     0     0     0     0   413   413
#       24      0     0     0     0     0     0     0   386     0     0     0     0     0     0   386
#       25      0     0     0     0     0     0     0     0     0     0   337     0     0     0   337
#       26      0     0     0     0     0     0     0   305     0     0     0     0     0     0   305
#       27      0     0     0     0     0     0     0   142     0     0     0     0     0     0   142
#       28      0     0     0     0     0     0     0     0     0     0     0    70     0     0    70
#       Sum 10770  8880  7652  6023  3704  2569  2278  1594  1060   850   766   625   497   413 47681

Sys.time()
g_k5 <- buildSNNGraph(sce, k = 5, use.dimred = 'PCA')
Sys.time()
## About 5 minutes

Sys.time()
g_walk_k5 <- igraph::cluster_walktrap(g_k5)
Sys.time()

## Takes about 33 min
# [1] "2019-11-14 14:35:36 EST"
# [1] "2019-11-14 15:03:00 EST"

clust_k5 <- sort_clusters(g_walk_k5$membership)

save(g_k5, g_walk_k5, file = 'rda_scran/g_k5.Rdata')
sce_image_grid(sce, clust_k5, 'pdf_scran/grid_SNN_k5_noXY.pdf', colors = cols)


## For it to be k = 7
clust_k5_k7 <- sort_clusters(igraph::cut_at(g_walk_k5, n = 7))
sce_image_grid(sce, clust_k5_k7, 'pdf_scran/grid_SNN_k5_k7_noXY.pdf')
addmargins(table(clust_k5_k7, sce$subject))
# clust_k5_k7 Br5292 Br5595 Br8100   Sum
#         1    13922   3069  10806 27797
#         2      427  11366    451 12244
#         3     2718    272    402  3392
#         4      520     24   1555  2099
#         5       84     38   1066  1188
#         6      342    504     77   923
#         7       20     11      7    38
#         Sum  18033  15284  14364 47681

addmargins(table(clust_k5_k7, clust_k50_k7))
#            clust_k50_k7
# clust_k5_k7     1     2     3     4     5     6     7   Sum
#         1   21690  3109  2686   105     0   114    93 27797
#         2     791 11347    92     1     0    12     1 12244
#         3      74    12  1126  2171     0     0     9  3392
#         4       0     0    13     0  1734   203   149  2099
#         5       0     5    73     1   176   933     0  1188
#         6       0     1   137     0     0     1   784   923
#         7       0     0    36     0     0     0     2    38
#         Sum 22555 14474  4163  2278  1910  1263  1038 47681




## From
## https://osca.bioconductor.org/clustering.html#assessing-cluster-separation

Sys.time()
ratio_k5 <- clusterModularity(g_k5, clust_k5, as.ratio = TRUE)
Sys.time()
ratio_k10 <- clusterModularity(g_k10, clust_k10, as.ratio = TRUE)
Sys.time()
ratio_k50 <- clusterModularity(g_k50, clust_k50, as.ratio = TRUE)
Sys.time()
save(ratio_k5, ratio_k10, ratio_k50, file = 'rda_scran/ratio_k5_10_50.Rdata')

pdf('pdf_scran/ratio_k5_10_50.pdf')
pheatmap(
    log2(ratio_k5 + 1),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("white", "blue"))(100)
)
pheatmap(
    log2(ratio_k10 + 1),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("white", "blue"))(100)
)
pheatmap(
    log2(ratio_k50 + 1),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("white", "blue"))(100)
)
dev.off()


cluster.gr_k5 <- igraph::graph_from_adjacency_matrix(ratio_k5,
    mode = "upper",
    weighted = TRUE,
    diag = FALSE)
Sys.time()
cluster.gr_k10 <- igraph::graph_from_adjacency_matrix(ratio_k10,
    mode = "upper",
    weighted = TRUE,
    diag = FALSE)
Sys.time()
cluster.gr_k50 <- igraph::graph_from_adjacency_matrix(ratio_k50,
    mode = "upper",
    weighted = TRUE,
    diag = FALSE)
Sys.time()
save(cluster.gr_k5, cluster.gr_k10, cluster.gr_k50, file = 'rda_scran/cluster.gr_k5_10_50.Rdata')


pdf('pdf_scran/cluster.gr_k5_10_50.pdf')
set.seed(11001010)
plot(cluster.gr_k5, edge.width = igraph::E(cluster.gr_k5)$weight * 1 / 4)
plot(cluster.gr_k10, edge.width = igraph::E(cluster.gr_k10)$weight * 1 /
        4)
plot(cluster.gr_k50, edge.width = igraph::E(cluster.gr_k50)$weight * 1 /
        4)
dev.off()


## Focus on the k50 SNN clusters
## From https://stackoverflow.com/questions/5812493/how-to-add-leading-zeros
clust_k50_d2 <-
    formatC(
        as.integer(clust_k50),
        width = 2,
        format = "d",
        flag = "0"
    )
cIndexes <- splitit(paste0(sce$sample_name, '_', clust_k50_d2))

## Adapted from collapse_clusters.R
## Collapse UMIs
umiComb <-
    sapply(cIndexes, function(ii)
        rowSums(assays(sce)$counts[top.hvgs, ii, drop = FALSE]))
dim(umiComb)
# [1] 1942  306

## Get a sample-specific size factors, instead of sample/cluster size factors
umiComb_sample <-
    sapply(splitit(ss(colnames(umiComb), '_', 1)), function(ii) {
        rowSums(umiComb[, ii, drop = FALSE])
    })
## Same as further below
# umiComb_sample_k50_k7 <- sapply(splitit(ss(colnames(umiComb_k50_k7), '_', 1)), function(ii) { rowSums(umiComb_k50_k7[, ii, drop = FALSE]) })
# identical(umiComb_sample, umiComb_sample_k50_k7)
umiComb_sample_size_fac <- librarySizeFactors(umiComb_sample)

umiComb_sample_size_fac_k50 <-
    rep(umiComb_sample_size_fac, lengths(splitit(ss(
        colnames(umiComb), '_', 1
    ))))
names(umiComb_sample_size_fac_k50) <- colnames(umiComb)

umiCombLog_k50 <-
    logNormCounts(SingleCellExperiment(list(counts = umiComb)), size_factors = umiComb_sample_size_fac_k50)
d_k50 <- dist(t(assays(umiCombLog_k50)$logcounts))
h_k50 <- hclust(d_k50)

pdf('pdf_scran/dendro_k50.pdf', width = 35)
palette(RColorBrewer::brewer.pal(12, 'Paired'))
myplclust(h_k50,
    labels = ss(names(cIndexes), '_', 2),
    lab.col = as.numeric(factor(ss(
        names(cIndexes), '_', 1
    ))))
dev.off()

cc_k50 <- cor(assays(umiCombLog_k50)$logcounts)
mean(cc_k50[upper.tri(cc_k50)])
# [1] 0.6951459

## Previously without the size factors
umiCombLog <-
    logNormCounts(SingleCellExperiment(list(counts = umiComb)))
d <- dist(t(assays(umiCombLog)$logcounts))
h <- hclust(d)

pdf('pdf_scran/dendro_k50_no_sizeFactors.pdf', width = 35)
palette(RColorBrewer::brewer.pal(12, 'Paired'))
myplclust(h,
    labels = ss(names(cIndexes), '_', 2),
    lab.col = as.numeric(factor(ss(
        names(cIndexes), '_', 1
    ))))
dev.off()

cc <- cor(assays(umiCombLog)$logcounts)
mean(cc[upper.tri(cc)])
# [1] 0.591435


## Build the annotation data.frame for the umi/cluster combination
col_df <- data.frame(cluster = factor(ss(names(cIndexes), '_', 2)),
    sample = ss(names(cIndexes), '_', 1))
col_df$subjpos <-
    sce$subject_position[match(as.character(col_df$sample), sce$sample_name)]
rownames(col_df) <- colnames(umiComb)

## Select colors myself
## Make it so the names have two digits
cols_d2 <- cols
names(cols_d2) <- unique(clust_k50_d2)

cols_sample <- RColorBrewer::brewer.pal(12, 'Paired')
names(cols_sample) <- unique(col_df$sample)

cols_subjpos <- RColorBrewer::brewer.pal(6, 'Paired')
names(cols_subjpos) <- unique(col_df$subjpos)

ann_colors <- list(cluster = cols_d2,
    sample = cols_sample,
    subjpos = cols_subjpos)

pdf('pdf_scran/pheatmap_umis_combined_k50_no_sizeFactors.pdf',
    height = 24)
pheatmap(
    assays(umiCombLog)$logcounts[top.hvgs,],
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("white", "blue"))(100),
    annotation_col = col_df,
    annotation_names_col = TRUE,
    annotation_colors = ann_colors,
    show_rownames = FALSE,
    show_colnames = FALSE
)
dev.off()


pdf('pdf_scran/pheatmap_umis_combined_k50.pdf', height = 24)
pheatmap(
    assays(umiCombLog_k50)$logcounts[top.hvgs,],
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("white", "blue"))(100),
    annotation_col = col_df,
    annotation_names_col = TRUE,
    annotation_colors = ann_colors,
    show_rownames = FALSE,
    show_colnames = FALSE
)
dev.off()



## Plot each cluster at a time in the grid, so we can clearly find where
## they are located
sce_image_grid_by_clus(sce,
    clust_k50_k7,
    'pdf_scran/grid_SNN_k50_k7_noXY_byCluster.pdf',
    ... = 'SNN k50 (cut to k7)')
sce_image_grid_by_clus(sce,
    clust_k50,
    'pdf_scran/grid_SNN_k50_noXY_byCluster.pdf',
    ... = 'SNN k50')



## Find marker genes https://osca.bioconductor.org/marker-gene-detection.html#using-pairwise-t-tests
Sys.time()
markers_wmw_k50 <-
    findMarkers(
        sce,
        clust_k50_d2,
        test = 'wilcox',
        block = sce$subject_position,
        direction = 'up'
    )
Sys.time()

## Takes about 11 minutes
# [1] "2019-11-18 14:51:48 EST"
# [1] "2019-11-18 15:02:18 EST"

save(markers_wmw_k50 , file = 'rda_scran/markers_wmw_k50.Rdata')


to_symbols <- function(x) {
    m <- match(rownames(x), rownames(rowData(sce)))
    rownames(x) <- rowData(sce)$gene_name[m]
    return(x)
}

AUCs_k50 <- lapply(markers_wmw_k50, function(interesting.wmw) {
    best.set <- interesting.wmw[interesting.wmw$Top <= 5, ]
    AUCs <- to_symbols(as.matrix(best.set[, -(1:3)]))
    colnames(AUCs) <- sub("AUC.", "", colnames(AUCs))
    return(AUCs)
})

## Adapted from
## https://osca.bioconductor.org/marker-detection.html#using-the-wilcoxon-rank-sum-test
## where they say:
## "A value greater than 0.5 indicates that the gene is upregulated in the current
## cluster compared to the other cluster, while values less than 0.5 correspond to
## downregulation. We would typically expect AUCs of 0.7-0.8 for a strongly
## upregulated candidate marker."
pdf('pdf_scran/AUCs_k50.pdf')
lapply(AUCs_k50, function(AUCs) {
    print(pheatmap(
        AUCs,
        breaks = seq(0, 1, length.out = 21),
        color = viridis::viridis(21)
    ))
    return(invisible(NULL))
})
dev.off()


Sys.time()
markers_binom_k50 <-
    findMarkers(
        sce,
        clust_k50_d2,
        test = 'binom',
        block = sce$subject_position,
        direction = 'up'
    )
Sys.time()
## Takes about 6 minutes
# [1] "2019-11-18 15:09:41 EST"
# [1] "2019-11-18 15:15:01 EST"

save(markers_binom_k50, file = 'rda_scran/markers_binom_k50.Rdata')


top_binom_k50 <-
    lapply(markers_binom_k50, function(binom) {
        head(rownames(binom), n = 6)
    })
sce$clust_k50_d2 <- clust_k50_d2


pdf('pdf_scran/top_binom_k50.pdf', width = 14)
lapply(top_binom_k50, function(topgenes) {
    p <- plotExpression(sce, x = 'clust_k50_d2', features = topgenes)
    ## Switch to symbols
    p$data$Feature <-
        rowData(sce)$gene_name[match(p$data$Feature, rowData(sce)$gene_id)]
    print(p)
    return(invisible(NULL))
})
dev.off()





#########################################################################
## Repeat some of the same as in the earlier part but with clust_k50_k7
#########################################################################


cIndexes_k50_k7 <-
    splitit(paste0(sce$sample_name, '_', clust_k50_k7))

## Adapted from collapse_clusters.R
## Collapse UMIs
umiComb_k50_k7 <-
    sapply(cIndexes_k50_k7 , function(ii)
        rowSums(assays(sce)$counts[top.hvgs, ii, drop = FALSE]))
dim(umiComb_k50_k7)
# [1] 1942  82

## With the sample-level size factors
umiComb_sample_size_fac_k50_k7 <-
    rep(umiComb_sample_size_fac, lengths(splitit(ss(
        colnames(umiComb_k50_k7), '_', 1
    ))))
names(umiComb_sample_size_fac_k50_k7) <- colnames(umiComb_k50_k7)

umiCombLog_k50_k7 <-
    logNormCounts(SingleCellExperiment(list(counts = umiComb_k50_k7)), size_factors = umiComb_sample_size_fac_k50_k7)

d_k50_k7 <- dist(t(assays(umiCombLog_k50_k7)$logcounts))
h_k50_k7 <- hclust(d_k50_k7)


pdf('pdf_scran/dendro_k50_k7.pdf', width = 30)
palette(RColorBrewer::brewer.pal(12, 'Paired'))
myplclust(h_k50_k7,
    labels = ss(names(cIndexes_k50_k7), '_', 2),
    lab.col = as.numeric(factor(ss(
        names(cIndexes_k50_k7), '_', 1
    ))))

palette(cols[1:7])
myplclust(h_k50_k7,
    labels = ss(names(cIndexes_k50_k7), '_', 1),
    lab.col = as.numeric(factor(ss(
        names(cIndexes_k50_k7), '_', 2
    ))))
dev.off()


cc_k50_k7 <- cor(assays(umiCombLog_k50_k7)$logcounts)
mean(cc_k50_k7[upper.tri(cc_k50_k7)])
# [1] 0.762597


## without the sample-level size factors
umiCombLog_k50_k7_noS <-
    logNormCounts(SingleCellExperiment(list(counts = umiComb_k50_k7)))

d_k50_k7_noS <- dist(t(assays(umiCombLog_k50_k7_noS)$logcounts))
h_k50_k7_noS <- hclust(d_k50_k7_noS)


pdf('pdf_scran/dendro_k50_k7_no_sizeFactors.pdf', width = 30)
palette(RColorBrewer::brewer.pal(12, 'Paired'))
myplclust(h_k50_k7_noS,
    labels = ss(names(cIndexes_k50_k7), '_', 2),
    lab.col = as.numeric(factor(ss(
        names(cIndexes_k50_k7), '_', 1
    ))))

palette(cols[1:7])
myplclust(h_k50_k7_noS,
    labels = ss(names(cIndexes_k50_k7), '_', 1),
    lab.col = as.numeric(factor(ss(
        names(cIndexes_k50_k7), '_', 2
    ))))
dev.off()

cc_k50_k7_noS <- cor(assays(umiCombLog_k50_k7_noS)$logcounts)
mean(cc_k50_k7_noS[upper.tri(cc_k50_k7_noS)])
# [1] 0.6735224




## Build the annotation data.frame for the umi/cluster combination
col_df_k50_k7 <- data.frame(cluster = factor(ss(names(cIndexes_k50_k7), '_', 2)),
    sample = ss(names(cIndexes_k50_k7), '_', 1))
col_df_k50_k7$subjpos <-
    sce$subject_position[match(as.character(col_df_k50_k7$sample), sce$sample_name)]
rownames(col_df_k50_k7) <- colnames(umiComb_k50_k7)


ann_colors_k50_k7 <- list(cluster = cols[seq_len(length(unique(clust_k50_k7)))],
    sample = cols_sample,
    subjpos = cols_subjpos)

pdf('pdf_scran/pheatmap_umis_combined_k50_k7.pdf', height = 24)
pheatmap(
    assays(umiCombLog_k50_k7)$logcounts[top.hvgs,],
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("white", "blue"))(100),
    annotation_col = col_df_k50_k7,
    annotation_names_col = TRUE,
    annotation_colors = ann_colors_k50_k7,
    show_rownames = FALSE,
    show_colnames = FALSE
)
dev.off()


pdf('pdf_scran/pheatmap_umis_combined_k50_k7_no_sizeFactors.pdf',
    height = 24)
pheatmap(
    assays(umiCombLog_k50_k7_noS)$logcounts[top.hvgs,],
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("white", "blue"))(100),
    annotation_col = col_df_k50_k7,
    annotation_names_col = TRUE,
    annotation_colors = ann_colors_k50_k7,
    show_rownames = FALSE,
    show_colnames = FALSE
)
dev.off()


## Find marker genes https://osca.bioconductor.org/marker-gene-detection.html#using-pairwise-t-tests
Sys.time()
markers_wmw_k50_k7 <-
    findMarkers(
        sce,
        clust_k50_k7,
        test = 'wilcox',
        block = sce$subject_position,
        direction = 'up'
    )
Sys.time()

## Takes about 7 minutes
# [1] "2019-11-20 14:19:45 EST"
# [1] "2019-11-20 14:26:36 EST"

save(markers_wmw_k50_k7 , file = 'rda_scran/markers_wmw_k50_k7.Rdata')


AUCs_k50_k7 <-
    lapply(markers_wmw_k50_k7, function(interesting.wmw) {
        best.set <- interesting.wmw[interesting.wmw$Top <= 5, ]
        AUCs <- to_symbols(as.matrix(best.set[, -(1:3)]))
        colnames(AUCs) <- sub("AUC.", "", colnames(AUCs))
        return(AUCs)
    })

pdf('pdf_scran/AUCs_k50_k7.pdf')
lapply(AUCs_k50_k7, function(AUCs) {
    print(pheatmap(
        AUCs,
        breaks = seq(0, 1, length.out = 21),
        color = viridis::viridis(21)
    ))
    return(invisible(NULL))
})
dev.off()


Sys.time()
markers_binom_k50_k7 <-
    findMarkers(
        sce,
        clust_k50_k7,
        test = 'binom',
        block = sce$subject_position,
        direction = 'up'
    )
Sys.time()
## Takes < 1 minute
# [1] "2019-11-20 16:25:53 EST"
# [1] "2019-11-20 16:26:35 EST"

save(markers_binom_k50_k7, file = 'rda_scran/markers_binom_k50_k7.Rdata')


top_binom_k50_k7 <-
    lapply(markers_binom_k50_k7, function(binom) {
        head(rownames(binom), n = 6)
    })
sce$clust_k50_k7 <- clust_k50_k7


pdf('pdf_scran/top_binom_k50_k7.pdf', width = 14)
lapply(top_binom_k50_k7, function(topgenes) {
    p <- plotExpression(sce, x = 'clust_k50_k7', features = topgenes)
    ## Switch to symbols
    p$data$Feature <-
        rowData(sce)$gene_name[match(p$data$Feature, rowData(sce)$gene_id)]
    print(p)
    return(invisible(NULL))
})
dev.off()







### For the SNN graph with K = 50, find which nested subset best matches
## the clusters from 10x Genomics labeled by Kristen Maynard and Keri Martinowich
clust_k5_list <- lapply(4:28, function(n) {
    message(paste(Sys.time(), 'n =', n))
    sort_clusters(igraph::cut_at(g_walk_k50, n = n))
})
names(clust_k5_list) <- paste0('k', 4:28)
save(clust_k5_list, file = 'rda_scran/clust_k5_list.Rdata')


library('readxl')

## Maynard
clust_10x_maynard_guess <- read_xlsx('guess_the_layer_Kristen.xlsx')
clust_10x_maynard_guess$sample_cluster <-
    with(clust_10x_maynard_guess, paste0(sample_name, '_', Cluster))

clust_10x <- paste0(sce$sample_name, '_', sce$Cluster)
clust_10x_layer_maynard <-
    clust_10x_maynard_guess$Layer[match(clust_10x, clust_10x_maynard_guess$sample_cluster)]
clust_10x_layer_maynard <-
    factor(clust_10x_layer_maynard,
        levels = c("1", "2_3", "4", "4_5", "5", "5_6", "6", "1_6", "WM"))

sce_image_grid_by_clus(
    sce,
    clust_10x_layer_maynard,
    'pdf_scran/grid_c10x_layer_maynard_byCluster.pdf',
    ... = 'c10x layers by Maynard'
)


cols_layers_maynard <-
    c(
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
    )[seq_len(length(unique(clust_10x_layer_maynard)))]
names(cols_layers_maynard) <- unique(clust_10x_layer_maynard)
sce_image_grid(
    sce,
    clust_10x_layer_maynard,
    'pdf_scran/grid_c10x_layer_maynard.pdf',
    sort_clust = FALSE,
    colors = cols_layers_maynard
)


## Keri
clust_10x_martinowich_guess <-
    read_xlsx('guess_the_layer_keri.xlsx')
clust_10x_martinowich_guess$sample_cluster <-
    with(clust_10x_martinowich_guess,
        paste0(sample_name, '_', Cluster))

clust_10x_layer_martinowich <-
    clust_10x_martinowich_guess$Layer[match(clust_10x, clust_10x_martinowich_guess$sample_cluster)]
clust_10x_layer_martinowich <-
    gsub('/', '_', clust_10x_layer_martinowich)
clust_10x_layer_martinowich <-
    factor(
        clust_10x_layer_martinowich,
        levels = c("1", "2", "2_3", "3", "4", "5", "5_6", "6", "1_6", "1_5", "WM")
    )

sce_image_grid_by_clus(
    sce,
    clust_10x_layer_martinowich,
    'pdf_scran/grid_c10x_layer_martinowich_byCluster.pdf',
    ... = 'c10x layers by Martinowich'
)


unique(clust_10x_layer_martinowich)[!unique(clust_10x_layer_martinowich) %in% unique(clust_10x_layer_maynard)]
# [1] 3   2   1_5
# Levels: 1 2 2_3 3 4 5 5_6 6 1_6 1_5 WM
## Re-use colors and add missing ones
cols_layers_martinowich <-
    c(cols_layers_maynard,
        '3' = "grey",
        '2' = "white",
        '1_5' = "purple")

sce_image_grid(
    sce,
    clust_10x_layer_martinowich,
    'pdf_scran/grid_c10x_layer_martinowich.pdf',
    sort_clust = FALSE,
    colors = cols_layers_martinowich
)



addmargins(
    table('Maynard' = clust_10x_layer_maynard,
        'Martinowich' = clust_10x_layer_martinowich)
)
#        Martinowich
# Maynard     1     2   2_3     3     4     5   5_6     6   1_6   1_5    WM   Sum
#     1    4665     0     0     0     0     0     0     0     0     0     0  4665
#     2_3  2706  1295  7829   954   238     0     0     0     0     0     0 13022
#     4       0     0     0   623  7059     0     0     0     0     0     0  7682
#     4_5     0     0     0     0     0   802     0     0     0     0     0   802
#     5       0     0     0     0  2102  2344     0   809     0     0     0  5255
#     5_6     0     0     0     0   469   363   752     0     0     0     0  1584
#     6       0     0     0     0     0  2813     0  4180     0     0     0  6993
#     1_6     0     0     0     0     0     0     0     0  1709   768     0  2477
#     WM      0     0     0     0     0     0     0     0     0     0  5201  5201
#     Sum  7371  1295  7829  1577  9868  6322   752  4989  1709   768  5201 47681
save(
    clust_10x_layer_maynard,
    clust_10x_layer_martinowich,
    cols_layers_maynard,
    cols_layers_martinowich,
    file = 'rda_scran/clust_10x_layer_maynard_martinowich.Rdata'
)

## Maynard
c10x_maynard_vs_snn_k50 <- lapply(clust_k5_list, function(cl) {
    table(clust_10x_layer_maynard , 'SNN k50' = cl)
})


addmargins(c10x_maynard_vs_snn_k50[[1]])
#                        SNN k50
# clust_10x_layer_maynard     1     2     3     4   Sum
#                     1    2565     4  2063    33  4665
#                     2_3 13014     1     7     0 13022
#                     4    7681     1     0     0  7682
#                     4_5   802     0     0     0   802
#                     5    5230    25     0     0  5255
#                     5_6  1570    13     1     0  1584
#                     6    6408   478     0   107  6993
#                     1_6  2314    35    77    51  2477
#                     WM   1608  2616   130   847  5201
#                     Sum 41192  3173  2278  1038 47681

pdf('pdf_scran/c10x_maynard_vs_snn_k50_rawSpotCounts.pdf')
lapply(c10x_maynard_vs_snn_k50, function(x) {
    rownames(x) <- paste(rownames(x), ';', rowSums(x))
    colnames(x) <- paste(colnames(x), ';', colSums(x))
    x[x == 0] <- NA
    print(pheatmap(
        x,
        color = viridis::viridis(21),
        cluster_rows = FALSE
    ))
    return(invisible(NULL))
})
dev.off()


pdf('pdf_scran/c10x_maynard_vs_snn_k50_fixedColorScale.pdf')
lapply(c10x_maynard_vs_snn_k50, function(x) {
    y <- sweep(x, 1, rowSums(x), '/') * 100
    rownames(y) <- paste(rownames(y), ';', rowSums(x))
    colnames(y) <- paste(colnames(y), ';', colSums(x))
    y[y == 0] <- NA
    print(pheatmap(
        y,
        breaks = seq(0, 100, length.out = 21),
        color = viridis::viridis(21),
        cluster_rows = FALSE
    ))
    return(invisible(NULL))
})
dev.off()


## Martinowich
c10x_martinowich_vs_snn_k50 <- lapply(clust_k5_list, function(cl) {
    table(clust_10x_layer_martinowich, 'SNN k50' = cl)
})


addmargins(c10x_martinowich_vs_snn_k50[[1]])
#                            SNN k50
# clust_10x_layer_martinowich     1     2     3     4   Sum
#                         1    5267     4  2067    33  7371
#                         2    1295     0     0     0  1295
#                         2_3  7826     1     2     0  7829
#                         3    1576     0     1     0  1577
#                         4    9865     3     0     0  9868
#                         5    6222    98     0     2  6322
#                         5_6   751     0     1     0   752
#                         6    4468   416     0   105  4989
#                         1_6  1579    28    56    46  1709
#                         1_5   735     7    21     5   768
#                         WM   1608  2616   130   847  5201
#                         Sum 41192  3173  2278  1038 47681

pdf('pdf_scran/c10x_martinowich_vs_snn_k50_rawSpotCounts.pdf')
lapply(c10x_martinowich_vs_snn_k50, function(x) {
    rownames(x) <- paste(rownames(x), ';', rowSums(x))
    colnames(x) <- paste(colnames(x), ';', colSums(x))
    x[x == 0] <- NA
    print(pheatmap(
        x,
        color = viridis::viridis(21),
        cluster_rows = FALSE
    ))
    return(invisible(NULL))
})
dev.off()


pdf('pdf_scran/c10x_martinowich_vs_snn_k50_fixedColorScale.pdf')
lapply(c10x_martinowich_vs_snn_k50, function(x) {
    y <- sweep(x, 1, rowSums(x), '/') * 100
    rownames(y) <- paste(rownames(y), ';', rowSums(x))
    colnames(y) <- paste(colnames(y), ';', colSums(x))
    y[y == 0] <- NA
    print(pheatmap(
        y,
        breaks = seq(0, 100, length.out = 21),
        color = viridis::viridis(21),
        cluster_rows = FALSE
    ))
    return(invisible(NULL))
})
dev.off()


## The second-round graphs:

pdf('pdf_scran/c10x_maynard_vs_snn_k50_rowPercent.pdf')
mapply(function(x, k) {
    y <- sweep(x, 1, rowSums(x), '/') * 100
    rownames(y) <- paste(rownames(y), ';', rowSums(x))
    colnames(y) <- paste(colnames(y), ';', colSums(x))
    y[y == 0] <- NA
    print(pheatmap(
        y,
        color = viridis::viridis(21),
        cluster_rows = FALSE,
        main = paste('Maynard vs SNN K50 cut at', k, '(row percents)')
    ))
    return(invisible(NULL))
},
    c10x_maynard_vs_snn_k50,
    names(c10x_maynard_vs_snn_k50))
dev.off()


pdf('pdf_scran/c10x_martinowich_vs_snn_k50_rowPercent.pdf')
mapply(function(x, k) {
    y <- sweep(x, 1, rowSums(x), '/') * 100
    rownames(y) <- paste(rownames(y), ';', rowSums(x))
    colnames(y) <- paste(colnames(y), ';', colSums(x))
    y[y == 0] <- NA
    print(pheatmap(
        y,
        color = viridis::viridis(21),
        cluster_rows = FALSE,
        main = paste('Martinowich vs SNN K50 cut at', k, '(row percents)')
    ))
    return(invisible(NULL))
},
    c10x_martinowich_vs_snn_k50,
    names(c10x_martinowich_vs_snn_k50))
dev.off()


pdf('pdf_scran/c10x_maynard_vs_snn_k50_colPercent.pdf')
mapply(function(x, k) {
    y <- sweep(x, 2, colSums(x), '/') * 100
    rownames(y) <- paste(rownames(y), ';', rowSums(x))
    colnames(y) <- paste(colnames(y), ';', colSums(x))
    y[y == 0] <- NA
    print(pheatmap(
        y,
        color = viridis::viridis(21),
        cluster_rows = FALSE,
        main = paste('Maynard vs SNN K50 cut at', k, '(col percents)')
    ))
    return(invisible(NULL))
},
    c10x_maynard_vs_snn_k50,
    names(c10x_maynard_vs_snn_k50))
dev.off()


pdf('pdf_scran/c10x_martinowich_vs_snn_k50_colPercent.pdf')
mapply(function(x, k) {
    y <- sweep(x, 2, colSums(x), '/') * 100
    rownames(y) <- paste(rownames(y), ';', rowSums(x))
    colnames(y) <- paste(colnames(y), ';', colSums(x))
    y[y == 0] <- NA
    print(pheatmap(
        y,
        color = viridis::viridis(21),
        cluster_rows = FALSE,
        main = paste('Martinowich vs SNN K50 cut at', k, '(col percents)')
    ))
    return(invisible(NULL))
},
    c10x_martinowich_vs_snn_k50,
    names(c10x_martinowich_vs_snn_k50))
dev.off()







## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

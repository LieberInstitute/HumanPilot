# screen -S sce
# qrsh -l mem_free=60G,h_vmem=60G,h_fsize=100G -pe local 4
# module load conda_R/3.6.x
library('SingleCellExperiment')
library('zinbwave')
library('clusterExperiment')
library('BiocParallel')
library('scran')
library('RColorBrewer')
library('sessioninfo')


dir.create('pdf_zinbwave', showWarnings = FALSE)
dir.create('rda_zinbwave', showWarnings = FALSE)


## From convert_sce.R
load('geom_spatial.Rdata', verbose = TRUE)

## From sce_scran.R
load('Human_DLPFC_Visium_processedData_sce_scran.Rdata', verbose = TRUE)


## For parallelization purposes
register(MulticoreParam(4))


table(rowSums(assay(sce) > 1))
table(rowSums(assay(sce) > 1) > 5)
# FALSE  TRUE
# 19352 14186

table(rowSums(assay(sce) > 2) > 5)
# FALSE  TRUE
# 24809  8729

table(rowSums(assay(sce) > 3) > 5)
# FALSE  TRUE
# 28660  4878

table(rowSums(assay(sce) > 5) > 5)
# FALSE  TRUE
# 31786  1752

# filter <- rowSums(assay(sce) > 1) > 5
# filtered <- sce[filter, ]


## makeFilterStats and zinbwave() doesn't work with the Matrix object
filtered <- sce[top.hvgs, ]
assays(filtered)$counts <- as.matrix(assays(filtered)$counts)
# assays(filtered)$logcounts <- as.matrix(assays(filtered)$logcounts)

## From
## https://www.bioconductor.org/packages/release/bioc/vignettes/scone/inst/doc/sconeTutorial.html#sample-filtering-with-metric_sample_filter
# > quantile(assay(sce)[assay(sce) > 0])
# 0%  25%  50%  75% 100%
#  1    1    1    2  610
num_reads <- quantile(assay(filtered)[assay(filtered) > 0])[4]
num_reads
# 75%
#   2
num_cells <- 0.25 * ncol(filtered)
num_cells
# [1] 11920.25
is_common <- rowSums(assay(filtered) >= num_reads ) >= num_cells
table(is_common)
# FALSE  TRUE
# 33039   499

## Drop in favor of the scran filtered list of genes output
# ## Continue with makeFilterStats()
# filtered <- makeFilterStats(filtered, filterStats="var", transFun = log1p)
# filtered <- filterData(filtered, percentile = 2000, filterStats="var")
# filtered
#
# table(rowSums(assay(filtered) > 1) > 5)
# #  TRUE
# # 2000
#
# table(rowSums(assay(filtered) >= num_reads ) >= num_cells)
# # FALSE  TRUE
# #  1501   499
 
## Adjust for sample
Sys.time()
clustered <- zinbwave(filtered, K = 50, X = "~ subject_position", residuals = TRUE,
    normalizedValues = TRUE, observationalWeights = TRUE, verbose = TRUE,
    BPPARAM = BiocParallel::MulticoreParam(4), epsilon = 1e12)
Sys.time()
## Takes about 20 hours to run!
# [1] "2019-11-12 14:11:12 EST"
# [1] "2019-11-13 09:44:34 EST"
save(clustered, file = 'rda_zinbwave/clustered.Rdata')
Sys.time()

## Set some colors
col_samples <- brewer.pal('Set3', n = length(unique(filtered$sample_name)))
names(col_samples) <- unique(filtered$sample_name)
    
## From
## https://bioconductor.github.io/BiocWorkshops/analysis-of-single-cell-rna-seq-data-dimensionality-reduction-clustering-and-lineage-inference.html#dimensionality-reduction
W <- reducedDim(clustered, 'zinbwave')
d <- dist(W)
length(d)
# [1] 1136715040

## Did not work: required too much mem
fit <- cmdscale(d, eig = TRUE, k = 2)

# pdf('pdf_zinbwave/mds_by_sample.pdf', useDingbats = FALSE)
# plot(fit$points, col = col_clus[filtered$sample_name], main = "",
#      pch = 20, xlab = "Component 1", ylab = "Component 2")
# legend(x = "topleft", legend = names(col_samples), cex = .5,
#     fill = col_samples, title = "Sample")
# dev.off()


## Try with scran
Sys.time()
g_k5 <- buildSNNGraph(clustered, k=5, use.dimred = 'zinbwave')
Sys.time()
## Takes about 2 minutes

Sys.time()
g_walk_k5 <- igraph::cluster_walktrap(g_k5)
Sys.time()

## Takes about 7 min
# [1] "2019-11-15 10:17:15 EST"
# [1] "2019-11-15 10:24:00 EST"

clust_k5 <- sort_clusters(g_walk_k5$membership)
length(unique(clust_k5))
# [1] 85
## Too many!

save(g_k5, g_walk_k5, file = 'rda_zinbwave/g_k5.Rdata')
# sce_image_grid(clustered, clust_k5, 'pdf_zinbwave/grid_SNN_k5_noXY.pdf', colors = cols)

## Try with K = 10
Sys.time()
g_k10 <- buildSNNGraph(clustered, k=10, use.dimred = 'zinbwave')
Sys.time()
## Takes about 2 minutes

Sys.time()
g_walk_k10 <- igraph::cluster_walktrap(g_k10)
Sys.time()
# [1] "2019-11-15 12:27:46 EST"
# [1] "2019-11-15 12:50:10 EST"

clust_k10 <- sort_clusters(g_walk_k10$membership)
length(unique(clust_k10))
# [1] 45

save(g_k10, g_walk_k10, file = 'rda_zinbwave/g_k10.Rdata')


## And with K = 50
Sys.time()
g_k50 <- buildSNNGraph(clustered, k=50, use.dimred = 'zinbwave')
Sys.time()
## Takes about 2 minutes

Sys.time()
g_walk_k50 <- igraph::cluster_walktrap(g_k50)
Sys.time()

clust_k50 <- sort_clusters(g_walk_k50$membership)
length(unique(clust_k50))

save(g_k50, g_walk_k50, file = 'rda_zinbwave/g_k50.Rdata')


## Remove since they are not needed right now
rm(d, W)

Sys.time()
## Fails due to memory
clustered <- RSEC(clustered, k0s = 4:15, alphas = c(0.1),
    betas = 0.8, reduceMethod="zinbwave",
    clusterFunction = "hierarchical01", minSizes=1,
    ncores = 1, isCount=FALSE,
    dendroReduce="zinbwave",
    subsampleArgs = list(resamp.num=100,
    clusterFunction="kmeans",
    clusterArgs=list(nstart=10)),
    verbose=TRUE,
    consensusProportion = 0.7,
    mergeMethod = "none", random.seed = 20191111,
    consensusMinSize = 10
)
Sys.time()


library('ClusterR')

km_rcpp <- function(x, k, checkArgs, cluster.only, ...) {
    km <- ClusterR::KMeans_rcpp(t(x), clusters = k, ...)
    if(cluster.only) {
        res <- km$clusters
    } else {
        res <- list(clustering = km$clusters)
    }
    return(res)
}
myCF <- ClusterFunction(clusterFUN=km_rcpp, inputType="X",algorithmType="K", outputType="vector")

clustered <- RSEC(clustered, k0s = 4:15, alphas = c(0.1),
    betas = 0.8, reduceMethod="zinbwave",
    clusterFunction = "hierarchical01", minSizes=1,
    ncores = 1, isCount=FALSE,
    dendroReduce="zinbwave",
    subsampleArgs = list(resamp.num=100,
    clusterFunction=myCF,
    clusterArgs=list(num_init=10)),
    verbose=TRUE,
    consensusProportion = 0.7,
    mergeMethod = "none", random.seed = 20191111,
    consensusMinSize = 10
)



## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

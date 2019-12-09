# Script to calculate clustering
# using small number of UMAP dimensions plus 2 spatial dimensions
# Lukas Weber, Dec 2019


library(SingleCellExperiment)
library(uwot)
library(scran)
library(scater)
library(ggplot2)
library(RColorBrewer)


# ---------
# load data
# ---------

# load scran output file (containing top 50 molecular PCs and 2 spatial coordinates)
load("../../data/Human_DLPFC_Visium_processedData_sce_scran.Rdata")

sce


# ------------
# extract data
# ------------

# extract PCs
dims_pcs <- reducedDim(sce, type = "PCA")

stopifnot(nrow(dims_pcs) == ncol(sce))

# extract spatial dimensions
dims_spatial <- colData(sce)[, c("imagerow", "imagecol")]

stopifnot(nrow(dims_spatial) == ncol(sce))


# -------------------------
# calculate UMAP dimensions
# -------------------------

# will aim to use top few (e.g. 5-10) UMAP dimensions (how much of the overall
# heterogeneity do these really capture?)

# calculate UMAP on PCs for faster runtime (could also calculate on all 1942
# highly variable genes instead for more accuracy)

# keep top 50 components
set.seed(123)
dims_umap <- umap(dims_pcs, n_components = 50)

stopifnot(nrow(dims_umap) == ncol(sce))


# ----------------
# scale dimensions
# ----------------

# need all dimensions (UMAP and spatial) to be on approximately comparable
# scales

# UMAP dimensions are already on a sensible scale, so can leave as is (note:
# don't do z-score scaling since this will scale up the less meaningful UMAP
# compenents)
summary(dims_umap)
mean(dims_umap[, 1])
sd(dims_umap[, 1])

range(dims_umap[, 1])
range(dims_umap[, 2])
range(dims_umap[, 3])

colnames(dims_umap) <- paste0("UMAP_", seq_len(ncol(dims_umap)))


# spatial dimensions: scale to e.g. min -5 and max 5, so they are on roughly
# similar scale as top few UMAP dimensions (note: z-score scaling doesn't really
# make sense for spatial coordinates)

# note: choice of these max and min values is very important! results will be
# highly sensitive to this
summary(as.data.frame(dims_spatial))
range(dims_spatial[, 1])
range(dims_spatial[, 2])

dims_spatial <- apply(as.matrix(dims_spatial), 2, function(col) {
    (col - min(col)) / (max(col) - min(col)) * 10 - 5
})

colnames(dims_spatial) <- c("spatial_x", "spatial_y")

summary(dims_spatial)

stopifnot(nrow(dims_spatial) == ncol(sce))


# ----------------------
# graph-based clustering
# ----------------------

# now can run standard Bioconductor graph-based clustering on subset of UMAP
# dimensions and scaled spatial dimensions

# number of UMAP dimensions to use
n_umap <- 5

dims_clus <- cbind(dims_umap[, seq_len(n_umap), drop = FALSE], dims_spatial)
head(dims_clus)


# see OSCA book (clustering chapter)
# note: use transpose
g <- buildSNNGraph(t(dims_clus), k = 10, d = ncol(dims_clus))
g_walk <- igraph::cluster_walktrap(g)

# default number of clusters (not using this for final results)
clus <- g_walk$membership
table(clus)

stopifnot(length(clus) == ncol(sce))

# choose number of clusters
clus <- igraph::cut_at(g_walk, n = 10)
table(clus)

stopifnot(length(clus) == ncol(sce))


# ------------
# plot results
# ------------

# display plot on original spatial coordinates

d_plot <- data.frame(
    x_coord = colData(sce)[, c("imagerow")], 
    y_coord = colData(sce)[, c("imagecol")], 
    cluster = as.factor(clus)
)

ggplot(d_plot, aes(x = x_coord, y = y_coord, color = cluster)) + 
    geom_point(size = 2, alpha = 0.5) + 
    coord_fixed() + 
    scale_color_brewer(palette = "Paired") + 
    theme_bw() + 
    ggtitle("Clustering on top few UMAP dims plus 2 spatial dims (scaled)")

ggsave("../plots/plot_clustering_umap_spatial.png", width = 7, height = 7)



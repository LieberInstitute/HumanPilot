###
library('SingleCellExperiment')
library('here')
library('jaffelab')
library('scater')
library('scran')
library('pheatmap')
library('readxl')
library('Polychrome')
library('cluster')
library('limma')
library('sessioninfo')
library('reshape2')

## Load data
load(here(
    'Analysis',
    'Human_DLPFC_Visium_processedData_sce_scran.Rdata'
))

## filter to variable genes
sce_hvg = sce[top.hvgs,]

## calc dist
dd = dist(t(assays(sce_hvg)$logcounts))
dd_mat = as.matrix(dd)

## too big to make long, do in chunks
sIndexes = split0(sce$sample_name)

## within slide distance
dist_within = t(sapply(sIndexes, function(ii) {
	cat(".")
	d_sub = dd_mat[ii,ii]
	c(mean = mean(d_sub), quantile(d_sub))
}))

## across slide distances
samplePairs = combn(length(sIndexes),2)

dist_across = apply(samplePairs, 2, function(i) {
	cat(".")
	ii = sIndexes[[i[1]]]
	jj = sIndexes[[i[2]]]
	d_sub = dd_mat[ii,jj]
	c(mean = mean(d_sub), quantile(d_sub))
})
save(distListPairs, file = "spotLevel_distances_allPairs.rda")


# save(dd, file = "distance_matrix_hvg_spotLevel.rda")
ddLong = melt(dd_mat)

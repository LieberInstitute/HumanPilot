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
library('lmerTest')

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
dist_across = t(dist_across)

###########
## combine and annotate
dist_within = as.data.frame(dist_within)
dist_within$Sample1 = rownames(dist_within)
dist_within$Sample2 = rownames(dist_within)

dist_across = as.data.frame(dist_across)
dist_across$Sample1 = names(sIndexes)[samplePairs[1,]]
dist_across$Sample2 = names(sIndexes)[samplePairs[2,]]

## add more metrics
metrics = colData(sce)
metrics = metrics[!duplicated(metrics$sample_name),]

dist_mat = rbind(dist_within, dist_across)

dist_mat$SubjPos1 = metrics$subject_position[match(dist_mat$Sample1, metrics$sample_name)]
dist_mat$SubjPos2 = metrics$subject_position[match(dist_mat$Sample2, metrics$sample_name)]
dist_mat$Subj1 = metrics$subject[match(dist_mat$Sample1, metrics$sample_name)]
dist_mat$Subj2 = metrics$subject[match(dist_mat$Sample2, metrics$sample_name)]

dist_mat$sameSlide = dist_mat$Sample1 == dist_mat$Sample2
dist_mat$spatialDup = dist_mat$SubjPos1 == dist_mat$SubjPos2
dist_mat$samePerson = dist_mat$Subj1 == dist_mat$Subj2
dist_mat$lab = factor(paste0(dist_mat$sameSlide, ":", dist_mat$spatialDup,":", dist_mat$samePerson))

## fix percentile names ##
colnames(dist_mat)[2:6] = paste0("perc_", seq(0,100,by=25))
dist_mat[dist_mat$lab == "TRUE:TRUE:TRUE",]

## mean distance
boxplot(mean ~ lab, data=dist_mat,outline=FALSE)
points(mean ~ jitter(as.numeric(lab),0.15), data=dist_mat,pch=21,bg="grey")

## median distance
boxplot(perc_50 ~ lab, data=dist_mat,outline=FALSE)
points(perc_50 ~ jitter(as.numeric(lab),0.15), data=dist_mat,pch=21,bg="grey")

summary(lm(perc_50 ~  samePerson + spatialDup + sameSlide , data=dist_mat))

# save(dd, file = "distance_matrix_hvg_spotLevel.rda")
ddLong = melt(dd_mat)

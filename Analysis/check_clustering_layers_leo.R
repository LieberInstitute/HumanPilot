## module load conda_R/3.6.x

## ----Libraries ------------------
library(parallel)
library(SummarizedExperiment)
library(Matrix)
library(RColorBrewer)
library(jaffelab) 
library(edgeR)

library('zinbwave')
library('SingleCellExperiment')
library('magrittr')
library('ggplot2')
library('Seurat')

## load rse list
load("Human_DLPFC_Visium_processedData_rseList.rda", verbose = TRUE)

## filter to expressed genes, lets be liberal
# exprsMat = sapply(rseList, function(rse) {
#     rowSums(assays(rse)$umis)
# })
#
# # hist(log2(rowSums(exprsMat)+1))
# table(rowSums(exprsMat) > 0)
# # liberal filter
# exprsIndex = which(rowSums(exprsMat) > 0)
# rseList = lapply(rseList, function(rse) rse[exprsIndex,])
#
# ## change rownames to symbols
# rseList = lapply(rseList, function(rse) {
#     rownames(rse) = rowData(rse)$gene_name
#     rse
# })


##
sample_n <- 9
sce <- SingleCellExperiment(
    assays = list(counts = assays(rseList[[sample_n]])$umis),
    rowData = rowData(rseList[[sample_n]]),
    colData = colData(rseList[[sample_n]]),
    metadata = metadata(rseList[[sample_n]])
)

## Relaxed filter from 5 to 1 from
## https://bioconductor.org/packages/release/bioc/vignettes/zinbwave/inst/doc/intro.html#gene-filtering
table(rowSums(assay(sce) >1))
table(rowSums(assay(sce) >1) > 5)
# FALSE  TRUE
# 25074  8464
filter <- rowSums(assay(sce) >1) > 5
sce <- sce[filter, ]

assay(sce) %>% as.matrix %>% log1p %>% rowVars -> vars
names(vars) <- rownames(sce)
vars <- sort(vars, decreasing = TRUE)
head(vars)


sce <- sce[names(vars)[seq_len(1e3)], ]

sce2 <- sce
## To avoid an error with t(dataY) at 
## https://github.com/drisso/zinbwave/blob/master/R/zinbwave.R#L254
assays(sce2)$counts <- as.matrix(assays(sce)$counts)
sce_zinb <- zinbwave(sce2, K = 2, epsilon=1000, verbose = TRUE,
    normalizedValues = TRUE, residuals = TRUE,
    observationalWeights = TRUE, BPPARAM = BiocParallel::MulticoreParam(4))

W <- reducedDim(sce_zinb)

data.frame(W, cluster=factor(colData(sce)$Cluster)) %>%
    ggplot(aes(W1, W2, colour=cluster)) + geom_point() + 
    scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()



## Seurat clustering
seu <- as.Seurat(x = sce_zinb, counts = "counts", data = "counts")
seu <- FindNeighbors(seu, reduction = "zinbwave",
    dims = 1:2 #this should match K,
)
seu <- FindClusters(object = seu, resolution = 0.2)

table(seu@meta.data$seurat_clusters)
#   1   2   3   4   5   6
# 924 843 795 776 531 357

identical(colnames(sce), as.character(seu@meta.data$barcode))

data.frame(W, cluster=factor(seu@meta.data$seurat_clusters)) %>%
    ggplot(aes(W1, W2, colour=cluster)) + geom_point() + 
    scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()
    
colData(sce)$seurat_clusters <- factor(seu@meta.data$seurat_clusters)
library(ggforce)
library(cowplot)
library(RColorBrewer)
library(grid)

# pdf("kmeans_spaceranger_grid.pdf",height=24, width=36)
# for(k in 2:10) {
#     cat(".")
    # plots_clusters = lapply(rseList, function(rse) {
        
        pdf('leo_sample9.pdf')
		d = as.data.frame(colData(sce))
        k <- 6        #
        # d$Group = d[,paste0("k", k)]
		ggplot(d, aes(x=imagecol,y=imagerow,fill=factor(seurat_clusters))) +
			geom_spatial(data=metadata(sce)$image,
				aes(grob=grob), x=0.5, y=0.5) +
			geom_point(shape = 21,  size = 1.25, stroke = 0.25)+
			coord_cartesian(expand=FALSE)+
			scale_fill_manual(values = c("#b2df8a","#e41a1c","#377eb8","#4daf4a","#ff7f00","gold",
				"#a65628", "#999999", "black", "grey", "white", "purple"))+
			xlim(0,max(sce$width)) +
			ylim(max(sce$height),0) +
			xlab("") + ylab("") +
			labs(fill = "Cluster")+
			guides(fill = guide_legend(override.aes = list(size=3)))+
			ggtitle(paste0(unique(sce$sample_name), " k=",k)) +
			theme_set(theme_bw(base_size = 10))+
			theme(panel.grid.major = element_blank(), 
					panel.grid.minor = element_blank(),
					panel.background = element_blank(), 
					axis.line = element_line(colour = "black"),
					axis.text = element_blank(),
					axis.ticks = element_blank())
        dev.off()
    # })

    # print(plot_grid(plotlist = plots_clusters))
# }
# dev.off()

#### k means
Ks = 6:12
names(Ks) = Ks
kList = lapply(Ks, function(k) {
	cat(k)
	mclapply(rseList, function(rse) {
		cat(".")
        # ## prep counts
        # u = assays(rse)$umis
        # lu = log2(u+1)
        # su = scale(lu)
        #
        # ## do kmeans on z-scale
        # km = kmeans(t(su), centers = k)
        
        hist(assays(sce_zinb)$normalizedValues[, 1], breaks = 50)
        hist(assays(sce_zinb)$residuals[, 1], breaks = 50)
        
        su <- assays(sce_zinb)$normalizedValues
        km <- kmeans(t(su), centers = k)
		
		## get DE by cluster
		cIndexes = splitit(km$cluster)
		tMat = sapply(cIndexes, function(ii) {
			x = rep(0,ncol(rse))
			x[ii] = 1
			f = lmFit(su, model.matrix(~x))
			eBayes(f)$t[,2]
		})
        cor(tMat)
        ## sample 9 output with k = 6
#             1          2          3           4          5          6
# 1  1.00000000 -0.8561794 -0.8809466  0.05557228 -0.7384702  0.4436690
# 2 -0.85617940  1.0000000  0.6637798  0.22706646  0.4649805 -0.5710342
# 3 -0.88094665  0.6637798  1.0000000 -0.44343799  0.9084135 -0.5829079
# 4  0.05557228  0.2270665 -0.4434380  1.00000000 -0.6351942  0.2017953
# 5 -0.73847020  0.4649805  0.9084135 -0.63519421  1.0000000 -0.5720947
# 6  0.44366904 -0.5710342 -0.5829079  0.20179525 -0.5720947  1.0000000

plot(tMat[, 1], tMat[, 2])
		
		list(kmeans = km, tMat = tMat)
	},mc.cores=4)
})


colData(sce)$k6_clusters <- factor(km$cluster)

pdf('leo_sample9_kmeans6.pdf')
d = as.data.frame(colData(sce))
k <- 6        #
# d$Group = d[,paste0("k", k)]
ggplot(d, aes(x=imagecol,y=imagerow,fill=factor(k6_clusters))) +
	geom_spatial(data=metadata(sce)$image,
		aes(grob=grob), x=0.5, y=0.5) +
	geom_point(shape = 21,  size = 1.25, stroke = 0.25)+
	coord_cartesian(expand=FALSE)+
	scale_fill_manual(values = c("#b2df8a","#e41a1c","#377eb8","#4daf4a","#ff7f00","gold",
		"#a65628", "#999999", "black", "grey", "white", "purple"))+
	xlim(0,max(sce$width)) +
	ylim(max(sce$height),0) +
	xlab("") + ylab("") +
	labs(fill = "Cluster")+
	guides(fill = guide_legend(override.aes = list(size=3)))+
	ggtitle(paste0(unique(sce$sample_name), " k=",k)) +
	theme_set(theme_bw(base_size = 10))+
	theme(panel.grid.major = element_blank(), 
			panel.grid.minor = element_blank(),
			panel.background = element_blank(), 
			axis.line = element_line(colour = "black"),
			axis.text = element_blank(),
			axis.ticks = element_blank())
dev.off()


sce$k6_clus1 <- as.integer(sce$k6_clusters == 1)
colData(sce2) <- colData(sce)

sce_zinb_k6_c1 <- zinbwave(sce2, K = 2, X = '~k6_clus1',
    epsilon=1000, verbose = TRUE,
    normalizedValues = TRUE, residuals = TRUE,
    observationalWeights = TRUE, BPPARAM = BiocParallel::MulticoreParam(4))
    
W_k6_c1 <- reducedDim(sce_zinb_k6_c1)

data.frame(W_k6_c1, cluster=factor(colData(sce)$k6_clusters)) %>%
    ggplot(aes(W1, W2, colour=cluster)) + geom_point() + 
    scale_color_brewer(type = "qual", palette = "Set1") + theme_classic()
    

hist(assays(sce_zinb_k6_c1)$normalizedValues[, 1], breaks = 50)


save(kList, file = "kmeans_k6to12_listOutput_withTstat.rda")
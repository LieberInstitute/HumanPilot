# qrsh -l mf=10G,h_vmem=11G,bluejay,h_fsize=100G -pe local 12
# cd /dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis
# module load conda_R/3.6.x

## ----Libraries ------------------
library(parallel)
library(SummarizedExperiment)
library(Matrix)
library(RColorBrewer)
library(jaffelab) 
library(edgeR)
library(flexclust) # for weighted kMeans

## load rse list
load("Human_DLPFC_Visium_processedData_rseList.rda")

## filter to expressed genes, lets be liberal
exprsMat = sapply(rseList, function(rse) {
	rowSums(assays(rse)$umis)
})

# hist(log2(rowSums(exprsMat)+1))
table(rowSums(exprsMat) > 0)
# liberal filter
exprsIndex = which(rowSums(exprsMat) > 0)
rseList = lapply(rseList, function(rse) rse[exprsIndex,])

## change rownames to symbols
rseList = lapply(rseList, function(rse) {
	rownames(rse) = rowData(rse)$gene_name
	rse
})
	
#### k means
Ks = 7:11
names(Ks) = Ks
kList = lapply(Ks, function(k) {
	cat(k)
	mclapply(rseList, function(rse) {
		cat(".")
		## prep counts
		u = assays(rse)$umis
		lu = log2(u+1)
		su = scale(lu,center=TRUE,scale=FALSE)

		## add x-y coordinates
		x = rbind(su, t(as.matrix(colData(rse)[,c("row", "col")])))
		wts = c(rep(1,nrow(u)), 10,10)
		
		## do kmeans on z-scale
		km = cclust(t(x), k = k, weights = wts, method="hardcl"	)
		
		## get DE by cluster
		cIndexes = splitit(clusters(km))
		tMat = sapply(cIndexes, function(ii) {
			x = rep(0,ncol(rse))
			x[ii] = 1
			f = lmFit(su, model.matrix(~x))
			f$coef[,2]
		})
		
		list(kmeans = km, tMat = tMat)
	},mc.cores=12)
})

save(kList, file = "kmeans_k7to11_listOutput_withTstat_weighted_10.rda")

################
##### plots ####
################

## add clusters to pheno data
for(i in seq(along=Ks)) {
	k = Ks[i]
	kk = kList[[i]]
	for(j in seq(along=rseList)) {
		rse = rseList[[j]]
		cl = clusters(kk[[j]]$kmeans)
		colData(rse) = cbind(colData(rse), cl)
		names(colData(rse))[ncol(colData(rse))] = paste0("k",k)
		rseList[[j]] = rse
	}
}

#### make plots
library(ggplot2)
library(ggforce)
library(cowplot)
pdf("kmeans_weighted_hardcl_grid_10.pdf",height=24, width=36)
for(k in 7:11) {
	cat(".")
	plots_clusters = lapply(rseList, function(rse) {
		d = as.data.frame(colData(rse))
		d$Group = d[,paste0("k", k)]
		ggplot(d, aes(x=imagecol,y=imagerow,fill=factor(Group))) +
			geom_spatial(data=metadata(rse)$image,
				aes(grob=grob), x=0.5, y=0.5) +
			geom_point(shape = 21,  size = 1.25, stroke = 0.25)+
			coord_cartesian(expand=FALSE)+
			scale_fill_manual(values = c("#b2df8a","#e41a1c","#377eb8","#4daf4a","#ff7f00","gold",
				"#a65628", "#999999", "black", "grey", "white", "purple"))+
			xlim(0,max(rse$width)) +
			ylim(max(rse$height),0) +
			xlab("") + ylab("") +
			labs(fill = "Cluster")+
			guides(fill = guide_legend(override.aes = list(size=3)))+
			ggtitle(paste0(unique(rse$sample_name), " k=",k)) +
			theme_set(theme_bw(base_size = 10))+
			theme(panel.grid.major = element_blank(), 
					panel.grid.minor = element_blank(),
					panel.background = element_blank(), 
					axis.line = element_line(colour = "black"),
					axis.text = element_blank(),
					axis.ticks = element_blank())
	})

	print(plot_grid(plotlist = plots_clusters))
}
dev.off()

###
## module load conda_R/3.6.x

## ----Libraries ------------------
## ----Libraries ------------------
library(tidyverse)
library(ggplot2)
library(Matrix)
library(Rmisc)
library(ggforce)
library(cowplot)
library(RColorBrewer)
library(grid)
library(SummarizedExperiment)
library(jaffelab)
library(parallel)

## load rse list
load("Human_DLPFC_Visium_processedData_rseList.rda")

## add kmeans
kmFiles= list.files("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X",
	pattern = "clustering", recur=TRUE, full=TRUE)
kmFiles = kmFiles[grep("kmeans", kmFiles)]

## get sample info
clusterDf = data.frame(SampleID = ss(kmFiles, "/", 8),
	K = as.numeric(ss(kmFiles, "_", 6)), 
	File = kmFiles, stringsAsFactors=FALSE)
clusterList = split(clusterDf, clusterDf$SampleID)
clusterList = clusterList[names(rseList)]

clList = lapply(clusterList, function(cc) {
	cc = cc[order(cc$K),]
	cl = lapply(cc$File, read.csv, 
		as.is=TRUE, row.names=1)
	clusters = do.call("cbind", cl)
	names(clusters) = paste0("k", cc$K)
	clusters
})

## join
rseList = mapply(function(rse, cl) {
	colData(rse) = cbind(colData(rse), cl)
	rse
}, rseList, clList)


## ---- fig.width = 16, fig.height = 8-------------------------------------

pdf("kmeans_spaceranger_grid.pdf",height=24, width=36)
for(k in 2:10) {
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

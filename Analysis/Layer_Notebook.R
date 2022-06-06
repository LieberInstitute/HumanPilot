## module load conda_R/3.6.x # devel


## ----Libraries ------------------
library(tidyverse)
library(ggplot2)
library(Matrix)
library(Rmisc)
library(ggforce)
library(rjson)
library(cowplot)
library(RColorBrewer)
library(grid)
library(readbitmap)
library(Seurat)
library(SummarizedExperiment)
library(rtracklayer)

## Function for plotting
geom_spatial <-  function(mapping = NULL,
                         data = NULL,
                         stat = "identity",
                         position = "identity",
                         na.rm = FALSE,
                         show.legend = NA,
                         inherit.aes = FALSE,
                         ...) {
  GeomCustom <- ggproto(
    "GeomCustom",
    Geom,
    setup_data = function(self, data, params) {
      data <- ggproto_parent(Geom, self)$setup_data(data, params)
      data
    },
    
    draw_group = function(data, panel_scales, coord) {
      vp <- grid::viewport(x=data$x, y=data$y)
      g <- grid::editGrob(data$grob[[1]], vp=vp)
      ggplot2:::ggname("geom_spatial", g)
    },
    
    required_aes = c("grob","x","y")
    
  )
  
  layer(
    geom = GeomCustom,
    mapping = mapping,
    data = data,
    stat = stat,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}



## get sample names
sample_names <- read.delim("lenas.txt", as.is=TRUE, header=FALSE)$V1
sample_names


##  10x output path
path = "/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/10X/"

## output
image_paths <- paste0(path, sample_names, "/tissue_lowres_image.png")
scalefactor_paths <- paste0(path, sample_names, "/scalefactors_json.json")
tissue_paths <- paste0(path, sample_names, "/tissue_positions_list.txt")
cluster_paths <- paste0(path, sample_names, "/", sample_names, "_analysis__clustering_graphclust_clusters.csv")
matrix_paths <- paste0(path, sample_names, "/", sample_names, "_filtered_feature_bc_matrix.h5")

all(file.exists(c(image_paths, scalefactor_paths, tissue_paths, cluster_paths, matrix_paths)))
# TRUE

## get annotation
map = read.delim("../10X/151675/151675_raw_feature_bc_matrix__features.tsv.gz",
	as.is=TRUE, header=FALSE,col.names=c("EnsemblID", "Symbol", "Type"))
## get GTF, this seems like what they used
gtf = import("/dcl01/ajaffe/data/lab/singleCell/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf")
gtf = gtf[gtf$type	== "gene"]
names(gtf) = gtf$gene_id
gtf = gtf[map$EnsemblID]
seqlevels(gtf)[1:25] = paste0("chr", seqlevels(gtf)[1:25])
mcols(gtf) = mcols(gtf)[,c(5:9)]

## ------------------------------------------------------------------------
images_cl <- lapply(image_paths, read.bitmap) 
dims = t(sapply(images_cl, dim))
colnames(dims) = c("height", "width", "channel")
dims = as.data.frame(dims)


## ------------------------------------------------------------------------
grobs <- lapply(images_cl, rasterGrob, width=unit(1,"npc"), height=unit(1,"npc"))
images_tibble <- tibble(sample=sample_names, grob=grobs)
images_tibble$height = dims$height
images_tibble$width = dims$width
images_tibble


## ------------------------------------------------------------------------
scales <- lapply(scalefactor_paths, function(x) fromJSON(file=x))


## ------------------------------------------------------------------------
clusters = lapply(cluster_paths, read.csv)
head(clusters[[1]])


## ------------------------------------------------------------------------
bcs <- list()
for (i in 1:length(sample_names)) {
   bcs[[i]] <- read.csv(tissue_paths[i],col.names=c("barcode","tissue","row","col","imagerow","imagecol"), header = FALSE)
   bcs[[i]]$sample_name <- sample_names[i]
   bcs[[i]]$imagerow <- bcs[[i]]$imagerow * scales[[i]]$tissue_lowres_scalef    # scale tissue coordinates for lowres image
   bcs[[i]]$imagecol <- bcs[[i]]$imagecol * scales[[i]]$tissue_lowres_scalef
   bcs[[i]]$tissue <- as.factor(bcs[[i]]$tissue)
   bcs[[i]] <- merge(bcs[[i]], clusters[[i]], by.x = "barcode", by.y = "Barcode", all = TRUE)
   bcs[[i]]$height <- images_tibble$height[i]
   bcs[[i]]$width <- images_tibble$width[i]
}

names(bcs) <- sample_names

head(bcs[[1]])


## ------------------------------------------------------------------------

## keep in regular genomics formatting
umiList <- lapply(matrix_paths, Read10X_h5)
names(umiList) = sample_names
sapply(umiList, dim)

rseList = mapply(function(u, bc) {
	rownames(bc) = bc$barcode
	bc = bc[colnames(u),c(1,7,2:6,8:ncol(bc))]
	rse = SummarizedExperiment(
		assays = list('umis' = u),
		rowRanges = gtf, colData = bc)
	rse$sum_umi = colSums(u)
	rse$sum_gene = colSums(u > 0)
	return(rse)
}, umiList, bcs)

## add images
for(i in seq(along=rseList)) {
	metadata(rseList[[i]])$image = images_tibble[i,]
}
## save out
save(rseList, geom_spatial,
	file = "Human_DLPFC_Visium_processedData_rseList.rda")


## ------------------------------------------------------------------------
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))


## ----  Just tissue per sample ------------------------------------

plots_tissue = lapply(rseList, function(rse) {
	d = as.data.frame(colData(rse))
	ggplot(d, aes(x=imagecol,y=imagerow,fill=sum_umi)) +
		geom_spatial(data=metadata(rse)$image,
			aes(grob=grob), x=0.5, y=0.5) +
		coord_cartesian(expand=FALSE)+
		scale_fill_gradientn(colours = myPalette(100))+
		xlim(0,max(rse$width)) +
		ylim(max(rse$height),0) +
		xlab("") + ylab("") +
		ggtitle(unique(rse$sample_name)) +
		theme_set(theme_bw(base_size = 10))+
		theme(panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_blank(), 
				axis.line = element_line(colour = "black"),
				axis.text = element_blank(),
				axis.ticks = element_blank())
})

pdf("example_tissue.pdf",height=24, width=36)
print(plot_grid(plotlist = plots_tissue))
dev.off()

#### by sample
pdf("tissue_by_sample.pdf")
for(i in seq(along=plots_tissue)) {
	print(plots_tissue[[i]])
}
dev.off()

## ----  Just tissue with blank spots per sample ------------------------------------
plots_tissueSpots = lapply(rseList, function(rse) {
	d = as.data.frame(colData(rse))
	ggplot(d, aes(x=imagecol,y=imagerow)) +
		geom_spatial(data=metadata(rse)$image,
			aes(grob=grob), x=0.5, y=0.5) +
		geom_point(shape = 21,  size = 1.25, stroke = 0.25)+
		coord_cartesian(expand=FALSE)+
		scale_fill_gradientn(colours = myPalette(100))+
		xlim(0,max(rse$width)) +
		ylim(max(rse$height),0) +
		xlab("") + ylab("") +
		ggtitle(unique(rse$sample_name)) +
		theme_set(theme_bw(base_size = 10))+
		theme(panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_blank(), 
				axis.line = element_line(colour = "black"),
				axis.text = element_blank(),
				axis.ticks = element_blank())
})

pdf("example_tissue_spotted.pdf",height=24, width=36)
print(plot_grid(plotlist = plots_tissueSpots))
dev.off()

#### by sample
pdf("tissueSpotted_by_sample.pdf")
for(i in seq(along=plots_tissueSpots)) {
	print(plots_tissueSpots[[i]])
}
dev.off()


## ----  UMIs per sample ------------------------------------
plots_umis = lapply(rseList, function(rse) {
	d = as.data.frame(colData(rse))
	ggplot(d, aes(x=imagecol,y=imagerow,fill=sum_umi)) +
		geom_spatial(data=metadata(rse)$image,
			aes(grob=grob), x=0.5, y=0.5) +
		geom_point(shape = 21,  size = 1.25, stroke = 0.25)+
		coord_cartesian(expand=FALSE)+
		scale_fill_gradientn(colours = myPalette(100))+
		xlim(0,max(rse$width)) +
		ylim(max(rse$height),0) +
		xlab("") + ylab("") +
		labs(fill = "Total UMI")+

		ggtitle(unique(rse$sample_name)) +
		theme_set(theme_bw(base_size = 10))+
		theme(panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_blank(), 
				axis.line = element_line(colour = "black"),
				axis.text = element_blank(),
				axis.ticks = element_blank())
})

pdf("example_umi.pdf",height=24, width=36)
print(plot_grid(plotlist = plots_umis))
dev.off()

#### by sample
pdf("umi_by_sample.pdf")
for(i in seq(along=plots_umis)) {
	print(plots_umis[[i]])
}
dev.off()


## ---- gene counts by sample ----------------------------------
plots_genes = lapply(rseList, function(rse) {
	d = as.data.frame(colData(rse))
	ggplot(d, aes(x=imagecol,y=imagerow,fill=sum_gene)) +
		geom_spatial(data=metadata(rse)$image,
			aes(grob=grob), x=0.5, y=0.5) +
		geom_point(shape = 21,  size = 1.25, stroke = 0.25)+
		coord_cartesian(expand=FALSE)+
		scale_fill_gradientn(colours = myPalette(100))+
		xlim(0,max(rse$width)) +
		ylim(max(rse$height),0) +
		xlab("") + ylab("") +
		labs(fill = "Total Gene")+
		ggtitle(unique(rse$sample_name)) +
		theme_set(theme_bw(base_size = 10))+
		theme(panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_blank(), 
				axis.line = element_line(colour = "black"),
				axis.text = element_blank(),
				axis.ticks = element_blank())
})


pdf("example_gene.pdf",height=24, width=36)
print(plot_grid(plotlist = plots_genes))
dev.off()

#### by sample
pdf("geneCount_by_sample.pdf")
for(i in seq(along=plots_genes)) {
	print(plots_genes[[i]])
}
dev.off()

## ---- fig.width = 16, fig.height = 8-------------------------------------
plots_clusters = lapply(rseList, function(rse) {
	d = as.data.frame(colData(rse))
	ggplot(d, aes(x=imagecol,y=imagerow,fill=factor(Cluster))) +
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
		ggtitle(unique(rse$sample_name)) +
		theme_set(theme_bw(base_size = 10))+
		theme(panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				panel.background = element_blank(), 
				axis.line = element_line(colour = "black"),
				axis.text = element_blank(),
				axis.ticks = element_blank())
})

pdf("example_cluster.pdf",height=24, width=36)
print(plot_grid(plotlist = plots_clusters))
dev.off()


#### by sample
pdf("cluster_by_sample.pdf")
for(i in seq(along=plots_clusters)) {
	print(plots_clusters[[i]])
}
dev.off()

#############################
## Layer and other markers
geneTab = read.csv("../mouse_layer_marker_info_cleaned.csv",as.is=TRUE)

library(biomaRt)
ensembl = useMart("ENSEMBL_MART_ENSEMBL", 
		dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
sym = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene"), 
	mart=ensembl)

geneTab$hgnc_symbol = sym$hgnc_symbol[match(geneTab$HumanEnsID, sym$ensembl_gene_id)]
geneTab$hgnc_symbol[geneTab$HumanEnsID == "ENSG00000275700"] = "AATF"

symbol = c("SLC17A7", "BDNF", "MBP", "MOBP", "GFAP", "MOG", "SNAP25", "GAD2", "CAMK2A", 
	"AQP4", "CD74", "FOXP2", "PDGFRA", "DLG4", geneTab$hgnc_symbol)
type = c("Excit", "Interest", "OLIGO", "OLIGO", "ASTRO", "OLIGO", "Neuron", "Inhib",
	"Neuron", "MICRO", "ASTRO", "NEURON", "OPC", "PSD95",paste("Layer", geneTab$Label))
dir.create("pdfs_grid")
dir.create("pdfs_single")

map = rowData(rseList[[1]])
ens = map$gene_id[match(symbol, map$gene_name)]


for(j in seq(along=symbol)) {
# for(j in 1:2) {
	cat(".")
	g = symbol[j]
	label = type[j]
	e = ens[j]
	
	plots = lapply(rseList, function(rse) {
		d = as.data.frame(colData(rse))
		d$UMI = assays(rse)$umis[e,]
		
		ggplot(d, aes(x=imagecol,y=imagerow,fill=UMI)) +
			geom_spatial(data=metadata(rse)$image,
				aes(grob=grob), x=0.5, y=0.5) +
			geom_point(shape = 21,  size = 1.25, stroke = 0.25)+
			coord_cartesian(expand=FALSE)+
			scale_fill_gradientn(colours = myPalette(100))+
			xlim(0,max(rse$width)) +
			ylim(max(rse$height),0) +
			xlab("") + ylab("") +
			labs(fill = "UMI")+
			ggtitle(paste(unique(rse$sample_name), g, label, sep=" - ")) +
			theme_set(theme_bw(base_size = 10))+
			theme(panel.grid.major = element_blank(), 
					panel.grid.minor = element_blank(),
					panel.background = element_blank(), 
					axis.line = element_line(colour = "black"),
					axis.text = element_blank(),
					axis.ticks = element_blank())
	})
	
	pdf(paste0("pdfs_grid/", g,"_", 
		gsub(" ", "", gsub("/","-", label)), ".pdf"),height=24, width=36)
	print(plot_grid(plotlist = plots))
	dev.off()
		
	#### by sample
	pdf(paste0("pdfs_single/", g,"_",
		gsub(" ", "", gsub("/","-", label)), "_bySample.pdf"))
	for(i in seq(along=plots)) {
		print(plots[[i]])
	}
	dev.off()	
}

## gene by cluster
# for(j in seq(along=symbol)) {
for(j in 1:2) {
	cat(".")
	g = symbol[j]
	label = type[j]
	e = ens[j]
	
	plots = lapply(rseList, function(rse) {
		d = as.data.frame(colData(rse))
		d$UMI = assays(rse)$umis[e,]
		d = d[d$UMI > 3,]
		
		ggplot(d, aes(x=imagecol,y=imagerow,fill=factor(Cluster))) +
			geom_spatial(data=metadata(rse)$image,
				aes(grob=grob), x=0.5, y=0.5) +
			geom_point(shape = 21,  size = 1.25, stroke = 0.25)+
			coord_cartesian(expand=FALSE)+
			scale_fill_manual(values = c("#b2df8a","#e41a1c","#377eb8","#4daf4a","#ff7f00","gold",
			"#a65628", "#999999", "black", "grey", "white", "purple"))+
			xlim(0,max(rse$width)) +
			ylim(max(rse$height),0) +
			xlab("") + ylab("") +
			labs(fill = "UMI")+
			ggtitle(paste(unique(rse$sample_name), g, label, sep=" - ")) +
			theme_set(theme_bw(base_size = 10))+
			theme(panel.grid.major = element_blank(), 
					panel.grid.minor = element_blank(),
					panel.background = element_blank(), 
					axis.line = element_line(colour = "black"),
					axis.text = element_blank(),
					axis.ticks = element_blank())
	})
	
	pdf(paste0("pdfs_grid/", g,"_", 
		gsub(" ", "", gsub("/","-", label)), "_sparse.pdf"),height=24, width=36)
	print(plot_grid(plotlist = plots))
	dev.off()
		
	#### by sample
	pdf(paste0("pdfs_single/", g,"_",
		gsub(" ", "", gsub("/","-", label)), "_bySample_sparse.pdf"))
	for(i in seq(along=plots)) {
		print(plots[[i]])
	}
	dev.off()	
}


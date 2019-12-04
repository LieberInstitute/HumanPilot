## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(
  cache = FALSE,
  cache.lazy = FALSE,
  tidy = TRUE
)


## ----Libraries, echo=TRUE, message=FALSE, warning=FALSE------------------
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


## ------------------------------------------------------------------------
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



## ----eval=FALSE, include=TRUE--------------------------------------------
## sample_names <- read.delim("lenas.txt", as.is=TRUE, header=FALSE)$V1
## sample_names


## ----eval=FALSE, include=TRUE--------------------------------------------
path = "/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X/"

## output
image_paths <- paste0(path, sample_names, "/tissue_lowres_image.png")
scalefactor_paths <- paste0(path, sample_names, "/scalefactors_json.json")
tissue_paths <- paste0(path, sample_names, "/tissue_positions_list.txt")
cluster_paths <- paste0(path, sample_names, "/", sample_names, "_analysis__clustering_graphclust_clusters.csv")
matrix_paths <- paste0(path, sample_names, "/", sample_names, "_filtered_feature_bc_matrix.h5")

all(file.exists(c(image_paths, scalefactor_paths, tissue_paths, cluster_paths, matrix_paths)))
# TRUE


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

matrix <- lapply(matrix_paths, Read10X_h5)
matrix = lapply(matrix, function(x) as.data.frame(t(x)))

head(matrix[[1]])


## ----message=FALSE, warning=FALSE----------------------------------------
umi_sum <- list() 

for (i in 1:length(sample_names)) {
  umi_sum[[i]] <- data.frame(barcode =  row.names(matrix[[i]]),
                             sum_umi = Matrix::rowSums(matrix[[i]]))
  
}
names(umi_sum) <- sample_names

umi_sum <- bind_rows(umi_sum, .id = "sample")
head(umi_sum)


## ----message=FALSE, warning=FALSE----------------------------------------
gene_sum <- list() 

for (i in 1:length(sample_names)) {
  gene_sum[[i]] <- data.frame(barcode =  row.names(matrix[[i]]),
                             sum_gene = Matrix::rowSums(matrix[[i]] != 0))
  
}
names(gene_sum) <- sample_names

gene_sum <- bind_rows(gene_sum, .id = "sample")
head(gene_sum)


## ------------------------------------------------------------------------
bcs_merge <- bind_rows(bcs, .id = "sample")
bcs_merge <- merge(bcs_merge,umi_sum, by = c("barcode", "sample"))
bcs_merge <- merge(bcs_merge,gene_sum, by = c("barcode", "sample"))
head(bcs_merge)


## ------------------------------------------------------------------------
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))


## ---- fig.width = 16, fig.height = 8-------------------------------------
plots <- list()

for (i in 1:length(sample_names)) {

plots[[i]] <- bcs_merge %>% 
  filter(sample ==sample_names[i]) %>% 
      ggplot(aes(x=imagecol,y=imagerow,fill=sum_umi)) +
                geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
                geom_point(shape = 21, colour = "black", size = 1.75, stroke = 0.5)+
                coord_cartesian(expand=FALSE)+
                scale_fill_gradientn(colours = myPalette(100))+
                xlim(0,max(bcs_merge %>% 
                            filter(sample ==sample_names[i]) %>% 
                            select(width)))+
                ylim(max(bcs_merge %>% 
                            filter(sample ==sample_names[i]) %>% 
                            select(height)),0)+
                xlab("") +
                ylab("") +
                ggtitle(sample_names[i])+
                labs(fill = "Total UMI")+
                theme_set(theme_bw(base_size = 10))+
                theme(panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(), 
                        axis.line = element_line(colour = "black"),
                        axis.text = element_blank(),
                        axis.ticks = element_blank())
}
pdf("example_umi.pdf",height=24, width=36)
print(plot_grid(plotlist = plots))
dev.off()


## ---- fig.width = 16, fig.height = 8-------------------------------------
plots <- list()

for (i in 1:length(sample_names)) {

plots[[i]] <- bcs_merge %>% 
  filter(sample ==sample_names[i]) %>% 
      ggplot(aes(x=imagecol,y=imagerow,fill=sum_gene)) +
                geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
                geom_point(shape = 21, colour = "black", size = 1.75, stroke = 0.5)+
                coord_cartesian(expand=FALSE)+
                scale_fill_gradientn(colours = myPalette(100))+
                xlim(0,max(bcs_merge %>% 
                            filter(sample ==sample_names[i]) %>% 
                            select(width)))+
                ylim(max(bcs_merge %>% 
                            filter(sample ==sample_names[i]) %>% 
                            select(height)),0)+
                xlab("") +
                ylab("") +
                ggtitle(sample_names[i])+
                labs(fill = "Total Genes")+
                theme_set(theme_bw(base_size = 10))+
                theme(panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(), 
                        axis.line = element_line(colour = "black"),
                        axis.text = element_blank(),
                        axis.ticks = element_blank())
}

pdf("example_gene.pdf",height=24, width=36)
print(plot_grid(plotlist = plots))
dev.off()


## ---- fig.width = 16, fig.height = 8-------------------------------------
plots <- list()

for (i in 1:length(sample_names)) {

plots[[i]] <- bcs_merge %>% 
  filter(sample ==sample_names[i]) %>%
  filter(tissue == "1") %>% 
      ggplot(aes(x=imagecol,y=imagerow,fill=factor(Cluster))) +
                geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
                geom_point(shape = 21, colour = "black", size = 1.75, stroke = 0.5)+
                coord_cartesian(expand=FALSE)+
                scale_fill_manual(values = c("#b2df8a","#e41a1c","#377eb8","#4daf4a","#ff7f00","gold", "#a65628", "#999999", "black", "grey", "white", "purple"))+
                xlim(0,max(bcs_merge %>% 
                            filter(sample ==sample_names[i]) %>% 
                            select(width)))+
                ylim(max(bcs_merge %>% 
                            filter(sample ==sample_names[i]) %>% 
                            select(height)),0)+
                xlab("") +
                ylab("") +
                ggtitle(sample_names[i])+
                labs(fill = "Cluster")+
                guides(fill = guide_legend(override.aes = list(size=3)))+
                theme_set(theme_bw(base_size = 10))+
                theme(panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(), 
                        axis.line = element_line(colour = "black"),
                        axis.text = element_blank(),
                        axis.ticks = element_blank())
}

pdf("example_cluster.pdf",height=24, width=36)
print(plot_grid(plotlist = plots))
dev.off()



## ---- fig.width = 16, fig.height = 8-------------------------------------
geneTab = read.csv("../mouse_layer_marker_info_cleaned.csv",as.is=TRUE)

library(biomaRt)
ensembl = useMart("ENSEMBL_MART_ENSEMBL", 
		dataset="hsapiens_gene_ensembl", host="feb2014.archive.ensembl.org")
sym = getBM(attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene"), 
	mart=ensembl)

geneTab$hgnc_symbol = sym$hgnc_symbol[match(geneTab$HumanEnsID, sym$ensembl_gene_id)]
geneTab$hgnc_symbol[geneTab$HumanEnsID == "ENSG00000275700"] = "AATF"

symbol = c("BDNF", "MBP", "MOBP", "GFAP", "MOG", "SNAP25", "GAD2", "CAMK2A", 
	"AQP4", "CD74", "FOXP2", "PDGFRA", "DLG4", geneTab$hgnc_symbol)
dir.create("pdfs")
for(j in seq(along=symbol)) {
	g = symbol[j]
	ge = enquo(g)

	plots <- list()

	for (i in 1:length(sample_names)) {

	plots[[i]] <- bcs_merge %>% 
					  filter(sample ==sample_names[i]) %>% 
					  bind_cols(select(matrix[[i]], g)) %>% 
	  ggplot(aes(x=imagecol,y=imagerow,fill=g)) +
					geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
					geom_point(shape = 21, colour = "black", size = 1.75, stroke = 0.5)+
					coord_cartesian(expand=FALSE)+
					scale_fill_gradientn(colours = myPalette(100))+
					xlim(0,max(bcs_merge %>% 
								filter(sample ==sample_names[i]) %>% 
								select(width)))+
					ylim(max(bcs_merge %>% 
								filter(sample ==sample_names[i]) %>% 
								select(height)),0)+
					xlab("") +
					ylab("") +
					ggtitle(paste(sample_names[i], g))+
					theme_set(theme_bw(base_size = 10))+
					theme(panel.grid.major = element_blank(), 
							panel.grid.minor = element_blank(),
							panel.background = element_blank(), 
							axis.line = element_line(colour = "black"),
							axis.text = element_blank(),
							axis.ticks = element_blank())
	}

	pdf(paste0("pdfs/", g, ".pdf"),height=24, width=36)
	print(plot_grid(plotlist = plots))
	dev.off()
}


## FEZF2
plots <- list()

for (i in 1:length(sample_names)) {

plots[[i]] <- bcs_merge %>% 
                  filter(sample ==sample_names[i]) %>% 
                  bind_cols(select(matrix[[i]], "FEZF2")) %>% 
  ggplot(aes(x=imagecol,y=imagerow,fill=FEZF2)) +
                geom_spatial(data=images_tibble[i,], aes(grob=grob), x=0.5, y=0.5)+
                geom_point(shape = 21, colour = "black", size = 1.75, stroke = 0.5)+
                coord_cartesian(expand=FALSE)+
                scale_fill_gradientn(colours = myPalette(100))+
                xlim(0,max(bcs_merge %>% 
                            filter(sample ==sample_names[i]) %>% 
                            select(width)))+
                ylim(max(bcs_merge %>% 
                            filter(sample ==sample_names[i]) %>% 
                            select(height)),0)+
                xlab("") +
                ylab("") +
                ggtitle(sample_names[i])+
                theme_set(theme_bw(base_size = 10))+
                theme(panel.grid.major = element_blank(), 
                        panel.grid.minor = element_blank(),
                        panel.background = element_blank(), 
                        axis.line = element_line(colour = "black"),
                        axis.text = element_blank(),
                        axis.ticks = element_blank())
}

pdf("example_FEZF2.pdf",height=24, width=36)
print(plot_grid(plotlist = plots))
dev.off()


## ---- out.width = "200px", echo=FALSE------------------------------------
knitr::include_graphics("~/public_html/Odin/Beta/example_notebook/hpca.jpg")


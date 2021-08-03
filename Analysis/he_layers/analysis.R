##
library(SummarizedExperiment)
library(limma)
library(recount)
library(jaffelab)
library(SingleCellExperiment)
library(here)
library(spatialLIBD)
library(RColorBrewer)
library(lattice)
library(pheatmap)

## load data
load("rse_gene_He_Layers_n102_annotated.Rdata")
## split by dataset
rse_gene_ds1 = rse_gene[,rse_gene$data_set == "DS1"]
rse_gene_ds2 = rse_gene[,rse_gene$data_set == "DS2"]
rse_gene_ds1$individuals = factor(rse_gene_ds1$individuals)

## reorder
rse_gene_ds1 = rse_gene_ds1[,order(rse_gene_ds1$individuals, rse_gene_ds1$section)]

## do pca
geneExprs_ds1 = log2(getRPKM(rse_gene_ds1,"Length")+1)
pca_ds1 = prcomp(t(geneExprs_ds1))
pcaVars_ds1 = getPcaVars(pca_ds1)

plot(pca_ds1$x, type="n",
	xlab = paste0("PC1: ", pcaVars_ds1[1], "% Var Expl"),
	ylab = paste0("PC2: ", pcaVars_ds1[1], "% Var Expl"))
text(pca_ds1$x[,1], pca_ds1$x[,2], rse_gene_ds1$section,
	col = as.numeric(factor(rse_gene$individuals)))

plot(pca_ds1$x[,1] ~ rse_gene_ds1$section,pch=21,
	bg=rse_gene_ds1$individuals, 
	ylab = paste0("PC1: ", pcaVars_ds1[1], "% Var Expl"))

## load our PCs
load("../Layer_Guesses/rda/sce_layer.Rdata", verbose = TRUE)
top.hvgs = top.hvgs[top.hvgs %in% rowData(rse_gene)$ensemblID]

pca = prcomp(t(assays(sce_layer)$logcounts[top.hvgs,]))
pcaVars =getPcaVars(pca)


## project their data
geneExprs_ds1_sub = geneExprs_ds1[match(top.hvgs,rowData(rse_gene)$ensemblID),]
geneExprs_ds1_scaled = scale(t(geneExprs_ds1_sub), pca$center, pca$scale)
genePCs_ds1_projected = geneExprs_ds1_scaled %*% pca$rotation 

plot(-1*pca$x[,1], pca$x[,2], pch = 21,cex=2,
	bg = libd_layer_colors[as.character(sce_layer$layer_guess)],
	xlab = paste0("PC1: ", pcaVars[1], "% Var Expl"),
	ylab = paste0("PC2: ", pcaVars[2], "% Var Expl"),
	ylim = c(-30,30), xlim = c(-65, 45))
palette(brewer.pal(6,"Dark2"))
text(-1*genePCs_ds1_projected[,1], genePCs_ds1_projected[,2],
	rse_gene_ds1$section,
	col = as.numeric(factor(rse_gene$individuals)))


#####################
## get our layer stats
######################
	
## load modeling outputs
load("../Layer_Guesses/rda/eb_contrasts.Rdata")
load("../Layer_Guesses/rda/eb0_list.Rdata")

## Extract the p-values
pvals0_contrasts <- sapply(eb0_list, function(x) {
    x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts) = rownames(eb_contrasts)
fdrs0_contrasts = apply(pvals0_contrasts, 2, p.adjust, "fdr")

## Extract the t-stats
t0_contrasts <- sapply(eb0_list, function(x) {
    x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts) = rownames(eb_contrasts)

#### line up to their expression data ###
mm = match(rownames(t0_contrasts), rowData(rse_gene_ds1)$ensemblID)

pvals0_contrasts = pvals0_contrasts[!is.na(mm),]
t0_contrasts = t0_contrasts[!is.na(mm),]
fdrs0_contrasts = fdrs0_contrasts[!is.na(mm),]
geneExprs_ds1_m = geneExprs_ds1[mm[!is.na(mm)],]

cor_t = cor(geneExprs_ds1_m,t0_contrasts)


### just layer specific genes from ones left
layer_specific_indices = mapply(function(t, p) {
    oo = order(t, decreasing = TRUE)[1:100]
},
    as.data.frame(t0_contrasts),
    as.data.frame(pvals0_contrasts))
layer_ind = unique(as.numeric(layer_specific_indices))

cor_t_layer = cor(geneExprs_ds1_m[layer_ind, ],
    t0_contrasts[layer_ind, ])
signif(cor_t_layer, 2)
cor_t_layer = cor_t_layer[,c(2:7,1)]

pdf("check_layer_corr.pdf",w=12,h=6)
par(mar=c(5,6,2,2),cex.axis=1.6,cex.lab=2,cex.main=1.8)
for(i in 1:ncol(cor_t_layer)){
	plot(rse_gene_ds1$section, cor_t_layer[,i], 
		ylim = c(-0.65,0.65),cex=2,
		ylab = paste0(colnames(cor_t_layer)[i]," (Correlation)"),
		xaxt = "n",xlab="Section",
		pch = 21, bg=factor(rse_gene_ds1$individuals))
	legend("top", levels(factor(rse_gene_ds1$individuals)), 
		pch =15, col=1:4,nc=4,cex=1.5)
	axis(1, at = 1:18)
}
dev.off()

## heatmap
cor_t_layer_wide =do.call("cbind", 
	split(as.data.frame(cor_t_layer), rse_gene_ds1$individuals))
cor_t_layer_wide_layer = cor_t_layer_wide[,rev(order(ss(colnames(cor_t_layer_wide), "\\.",2)))]
cor_t_layer_wide_layer = as.matrix(cor_t_layer_wide_layer)
rownames(cor_t_layer_wide_layer) = rse_gene_ds1$section[match(rownames(cor_t_layer_wide_layer), colnames(rse_gene_ds1))]

cor_t_layer_wide_person = cor_t_layer_wide[,rev(order(ss(colnames(cor_t_layer_wide), "\\.")))]
cor_t_layer_wide_person = as.matrix(cor_t_layer_wide_person)
rownames(cor_t_layer_wide_person) = rse_gene_ds1$section[match(rownames(cor_t_layer_wide_person), colnames(rse_gene_ds1))]

### heatmap
theSeq = seq(-.85, .85, by = 0.01)
my.col <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq))

pdf("he_sections_overlap_heatmap.pdf", width = 8)
print(
    levelplot(
        cor_t_layer_wide_layer,
        aspect = "fill",
        at = theSeq,
        col.regions = my.col,
        ylab = "",
        xlab = "Section",
        scales = list(x = list(cex = 1.5), y = list(cex = 1.5))
    )
)
print(
    levelplot(
        cor_t_layer_wide_layer_clust,
        aspect = "fill",
        at = theSeq,
        col.regions = my.col,
        ylab = "",
        xlab = "Section",
        scales = list(x = list(cex = 1.5), y = list(cex = 1.5))
    )
)

print(
    levelplot(
        cor_t_layer_wide_person,
        aspect = "fill",
        at = theSeq,
        col.regions = my.col,
        ylab = "",
        xlab = "Section",
        scales = list(x = list(cex = 1.5), y = list(cex = 1.5))
    )
)


dev.off()

anno = data.frame(Subj = ss(colnames(cor_t_layer_wide_person), "\\."),
	Layer = gsub("ayer", "", ss(colnames(cor_t_layer_wide_person), "\\.",2)),
	stringsAsFactors=FALSE)
anno$Layer = factor(anno$Layer, levels = c(paste0("L", 1:6), "WM"))
rownames(anno) = colnames(cor_t_layer_wide_person)
ann_colors = list(Layer = libd_layer_colors[1:7])
names(ann_colors$Layer ) = gsub("ayer", "",names(ann_colors$Layer ))

pdf("he_sections_overlap_pheatmap.pdf", width = 8)
pheatmap(cor_t_layer_wide_person, color = my.col,
	breaks = theSeq, fontsize = 16,
	annotation_col = anno, annotation_colors = ann_colors)
dev.off()

## make long version
pdf("he_sections_overlap_pheatmap_long.pdf", width = 14,h=7)
row_anno = as.data.frame(colData(rse_gene_ds1)[,c("section", "individuals")])
pheatmap(t(cor_t_layer), color = my.col,
	breaks = theSeq, fontsize = 16,
	annotation_col = row_anno)
dev.off()

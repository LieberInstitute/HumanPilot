##
library(SummarizedExperiment)
library(recount)
library(jaffelab)
library(edgeR)
library(SingleCellExperiment)
library(spatialLIBD)

## load counts
load("rse_exon_layerLevel_n76.Rdata")
load("rse_jx_layerLevel_n76.Rdata")
load("rse_gene_layerLevel_n76.Rdata")

## other phenotype data
load("Layer_Guesses/rda/sce_layer.Rdata")
colnames(sce_layer) = gsub("ayer", "", colnames(sce_layer))
sce_layer = sce_layer[,colnames(rse_gene)]
## add pheno
rse_gene$subject = sce_layer$subject
rse_gene$subject_position = sce_layer$subject_position

boxplot(totalAssignedGene ~ Layer, data = colData(rse_gene))
boxplot(overallMapping ~ Layer, data = colData(rse_gene))

#### check same pca as before ####
geneRpkm = getRPKM(rse_gene, "Length")
geneExprs = log2(geneRpkm+1)
gIndex = rowMeans(geneRpkm) > 0.1

## do pca
pca = prcomp(t(geneExprs[gIndex,]))
pcaVars = getPcaVars(pca)
pcaVars[1:10]

palette(libd_layer_colors[1:7])
plot(pca$x, pch =21, bg=factor(rse_gene$Layer),cex=2,
	xlab= paste0("PC1: ", pcaVars[1], "% Var Expl"),
	ylab= paste0("PC2: ", pcaVars[2], "% Var Expl"))
legend("bottomleft", levels(factor(rse_gene$Layer)),
	pch = 15, col = 1:7)

plot(pca$x[,1] ~ rse_gene$totalAssignedGene)
	
###################
### GENE LEVELS ###
###################

## voom
dge = DGEList(counts = assays(rse_gene[gIndex,])$counts, 
	genes = rowData(rse_gene[gIndex,]))
dge = calcNormFactors(dge)

## mean-variance
mod = model.matrix(~Layer + totalAssignedGene,
	data=colData(rse_gene))
vGene = voom(dge,mod,plot=TRUE)

gene_dupCorr = duplicateCorrelation(vGene$E, mod, 
	block=colData(rse_gene)$subject_position)
	
## Next for each layer test that layer vs the rest
layer_idx <- splitit(rse_gene$Layer)

eb0_list <- lapply(layer_idx, function(x) {
    res <- rep(0, ncol(rse_gene))
    res[x] <- 1
    m <- model.matrix(~ res + totalAssignedGene, data = colData(rse_gene))
    eBayes(
        lmFit(
            vGene$E,
            design = m,
            block = rse_gene$subject_position,
            correlation = gene_dupCorr$consensus.correlation
        )
    )
})

## Extract the p-values
pvals0_contrasts <- sapply(eb0_list, function(x) {
    x$p.value[, 2, drop = FALSE]
})
fdrs0_contrasts = apply(pvals0_contrasts, 2, p.adjust, "fdr")

t0_contrasts <- sapply(eb0_list, function(x) {
    x$t[, 2, drop = FALSE]
})

data.frame(
    'FDRsig' = colSums(apply(pvals0_contrasts, 2, p.adjust, 'fdr') < 0.05),
    'Pval10-6sig' = colSums(pvals0_contrasts < 1e-6),
    'Pval10-8sig' = colSums(pvals0_contrasts < 1e-8)
)

geneSigList = mapply(function(q, t) {
	rownames(rse_gene[gIndex,])[which(q < 0.05 & t > 0)]
}, as.data.frame(fdrs0_contrasts), as.data.frame(t0_contrasts))
lengths(geneSigList)
  # L1   L2   L3   L4   L5   L6   WM
# 1212  192   83   77  120  196 1691

#### compare to exons ?


###################
### exon LEVELS ###
###################

## reorder
rse_exon = rse_exon[order(rse_exon),]

## means by layer
exonRpkm = getRPKM(rse_exon, "Length")
exonLogMeans = sapply(layer_idx, function(ii) rowMeans(log2(exonRpkm[,ii] + 1)))

rowData(rse_exon) = cbind(rowData(rse_exon) , exonLogMeans)
rowData(rse_exon)$Overall = rowMeans(log2(exonRpkm+1))

## split by gene
emap = rowRanges(rse_exon)
emap_bygene = GRangesList(split(emap, factor(emap$ensemblID,
	levels = rownames(rse_gene))))
emap_bygene_filter = emap_bygene[rownames(rse_gene)[gIndex]]

## drop super long exons
emap_bygene_filter = emap_bygene_filter[width(emap_bygene_filter) < 2000]

## first vs last
last_vs_first_frac = sapply(emap_bygene_filter, function(x) {
	r = x$Overall[length(x)]/x$Overall[1] 
	if(all(strand(x) == "-")) r = 1/r
	return(r)
})
last_vs_first_frac = as.numeric(last_vs_first_frac)
last_vs_first_frac[!is.finite(last_vs_first_frac)] = NA
hist(log2(last_vs_first_frac))
num_exons = lengths(emap_bygene_filter)
gene_mean = log2(rowMeans(geneRpkm)[names(emap_bygene_filter)]+1)

plot(gene_mean, log2(last_vs_first_frac))

###### interpolate mean to 100 units

emap_bygene_filter_len = emap_bygene_filter[lengths(emap_bygene_filter) > 1]
cover = t(sapply(emap_bygene_filter_len, function(x) {
	x = x[order(x)] # put in order
	z= approx(as.numeric(x$Overall), n = 20)$y
	if(all(strand(x) == "-")) z = rev(z)
	return(z)
}))

## 
emap_bygene[rownames(rse_gene)[rowData(rse_gene)$Symbol == "BDNF"]]

cover_prop = prop.table(cover,1)
boxplot(cover_prop)

rse_exon_list = split(

eIndex = rowMeans(getRPKM(rse_exon, "Length")) > 0.1

## voom
dee = DGEList(counts = assays(rse_exon[eIndex,])$counts, 
	genes = rowData(rse_exon[eIndex,]))
dee = calcNormFactors(dee)

vExon = voom(dee,mod,plot=TRUE)
exon_dupCorr = duplicateCorrelation(vExon$E, mod, 
	block=colData(rse_gene)$subject_position)
	
eb0_list_exon <- lapply(layer_idx, function(x) {
    res <- rep(0, ncol(rse_gene))
    res[x] <- 1
    m <- model.matrix(~ res + totalAssignedGene, data = colData(rse_gene))
    eBayes(
        lmFit(
            vExon$E,
            design = m,
            block = rse_gene$subject_position,
            correlation = exon_dupCorr$consensus.correlation
        )
    )
})

## Extract the p-values
pvals0_contrasts_exons <- sapply(eb0_list_exon, function(x) {
    x$p.value[, 2, drop = FALSE]
})
fdrs0_contrasts_exons = apply(pvals0_contrasts_exons, 2, p.adjust, "fdr")

t0_contrasts_exons <- sapply(eb0_list_exon, function(x) {
    x$t[, 2, drop = FALSE]
})

## checks
data.frame(
    'FDRsig' = colSums(apply(pvals0_contrasts_exons, 2, p.adjust, 'fdr') < 0.05),
    'Pval10-6sig' = colSums(pvals0_contrasts_exons < 1e-6),
    'Pval10-8sig' = colSums(pvals0_contrasts_exons < 1e-8)
)

exonSigList = mapply(function(q, t) {
	unique(rowData(rse_exon[eIndex,])$ensemblID[which(q < 0.05 & t > 0)])
}, as.data.frame(fdrs0_contrasts_exons), as.data.frame(t0_contrasts_exons))
lengths(exonSigList)
  # L1   L2   L3   L4   L5   L6   WM
 # 974  179 1188   63  105  109 1930


mapply(function(e,g) {
	table(e %in% g)
}, exonSigList, geneSigList)
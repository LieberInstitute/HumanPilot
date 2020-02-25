###
#  module load conda_R/3.6.x
library(jaffelab)
library(Seurat)
library(scater)
library(DropletUtils)
library(limma)
library(lattice)
library(RColorBrewer)
library(pheatmap)

## read in sce.dlpfca
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_DLPFC-n2_cleaned-combined_SCE_MNTFeb2020.rda")

## drop ambig clusters
sce.dlpfc = sce.dlpfc[,sce.dlpfc$cellType != "Ambig.lowNtrxts"]

## numbers for paper
dim(sce.dlpfc)
length(unique(sce.dlpfc$prelimCluster))

## get pseudobulk
sce.dlpfc$PseudoSample = paste0(sce.dlpfc$sample,
			":", sce.dlpfc$prelimCluster)

## sum counts
cIndexes = splitit(sce.dlpfc$PseudoSample)
umiComb <- sapply(cIndexes, function(ii)
    rowSums(assays(sce.dlpfc)$counts[, ii, drop = FALSE]))

## filter pheno
phenoComb = colData(sce.dlpfc)[!duplicated(sce.dlpfc$PseudoSample), 
	c("prelimCluster", "collapsedCluster", "cellType", "PseudoSample")]
rownames(phenoComb) = phenoComb$PseudoSample
phenoComb = phenoComb[colnames(umiComb), ]
phenoComb = DataFrame(phenoComb)
phenoComb$prelimCluster = droplevels(phenoComb$prelimCluster)

sce_pseudobulk <-
    logNormCounts(SingleCellExperiment(
        list(counts = umiComb),
        colData = phenoComb,
        rowData = rowData(sce.dlpfc)
    ))
save(sce_pseudobulk, file = "rda/dlpfc_snRNAseq_pseudobulked.Rdata")

###############################
## extract expression
load("rda/dlpfc_snRNAseq_pseudobulked.Rdata")

mat <- assays(sce_pseudobulk)$logcounts

## Build a group model
mod <- with(colData(sce_pseudobulk),
    model.matrix(~ 0 + prelimCluster))
colnames(mod) <- gsub('prelimCluster', '', colnames(mod))

## get duplicate correlation
corfit <- duplicateCorrelation(mat, mod,
    block = sce_pseudobulk$sample)
save(corfit, file = "rda/dlpfc_snRNAseq_pseudobulked_dupCor.Rdata")

## Next for each layer test that layer vs the rest
cell_idx <- splitit(sce_pseudobulk$prelimCluster)

eb0_list_cell <- lapply(cell_idx, function(x) {
    res <- rep(0, ncol(sce_pseudobulk))
    res[x] <- 1
    m <- with(colData(sce_pseudobulk),
        model.matrix(~ res))
    eBayes(
        lmFit(
            mat,design = m,
            block = sce_pseudobulk$sample,
            correlation = corfit$consensus.correlation
        )
    )
})
save(eb0_list_cell, file = "rda/dlpfc_snRNAseq_pseudobulked_specific_Ts.Rdata")

##########
## Extract the p-values
load("rda/dlpfc_snRNAseq_pseudobulked_specific_Ts.Rdata")

pvals0_contrasts_cell <- sapply(eb0_list_cell, function(x) {
    x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts_cell) = rownames(mat)

t0_contrasts_cell <- sapply(eb0_list_cell, function(x) {
    x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts_cell) = rownames(mat)
fdrs0_contrasts_cell = apply(pvals0_contrasts_cell, 2, p.adjust, 'fdr')

data.frame(
    'FDRsig' = colSums(fdrs0_contrasts_cell < 0.05 &
            t0_contrasts_cell > 0),
    'Pval10-6sig' = colSums(pvals0_contrasts_cell < 1e-6 &
            t0_contrasts_cell > 0),
    'Pval10-8sig' = colSums(pvals0_contrasts_cell < 1e-8 &
            t0_contrasts_cell > 0)
)

   # FDRsig Pval10.6sig Pval10.8sig
# 1    1685         825         600
# 2     330         143         103
# 3     313         152         126
# 4     291         107          64
# 5     855         292         163
# 6     241          83          42
# 7     949         282         160
# 8     311         184         142
# 9     305         131          94
# 10    377         142          94
# 11    211          86          46
# 12    349         158         121
# 13    301         173         130
# 14    267         130         100
# 15    190          91          42
# 16    288         104          61
# 17    248         178          31
# 18    185          91          66
# 19    394         140          88
# 20    195          51          14
# 21    169          51          22
# 22    222          75          32
# 23    183          68          36
# 25    384         147          75
# 26    230          67          22
# 27    251         101          60
# 28    267         116          82
# 29    173          60          32
# 30    303         103          53
# 31    380         147          94


############################
### correlate to layer?? ###
############################

###################
## load modeling outputs
load("rda/eb_contrasts.Rdata")
load("rda/eb0_list.Rdata")

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

############
# line up ##

mm = match(rownames(pvals0_contrasts), rowData(sce_pseudobulk)$ID)

pvals0_contrasts = pvals0_contrasts[!is.na(mm), ]
t0_contrasts = t0_contrasts[!is.na(mm), ]
fdrs0_contrasts = fdrs0_contrasts[!is.na(mm), ]

pvals0_contrasts_cell = pvals0_contrasts_cell[mm[!is.na(mm)], ]
t0_contrasts_cell = t0_contrasts_cell[mm[!is.na(mm)], ]
fdrs0_contrasts_cell = fdrs0_contrasts_cell[mm[!is.na(mm)], ]

cor_t = cor(t0_contrasts_cell, t0_contrasts)
signif(cor_t, 2)

### just layer specific genes from ones left
layer_specific_indices = mapply(function(t, p) {
    oo = order(t, decreasing = TRUE)[1:100]
},
    as.data.frame(t0_contrasts),
    as.data.frame(pvals0_contrasts))
layer_ind = unique(as.numeric(layer_specific_indices))

cor_t_layer = cor(t0_contrasts_cell[layer_ind, ],
    t0_contrasts[layer_ind, ])
signif(cor_t_layer, 3)

### heatmap
theSeq = seq(-.85, .85, by = 0.01)
my.col <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq))

ct = colData(sce_pseudobulk)
ct = ct[!duplicated(sce_pseudobulk$prelimCluster),]
ct = ct[order(ct$cellType, ct$prelimCluster),]
ct$cellType = as.character(ct$cellType)
ct$cellType[ct$cellType == "Ambig.lowNtrxts"] ="Drop"
ct$lab = paste0(ct$prelimCluster, " (", ct$cellType,")")


dd = dist(1-cor_t_layer)
hc = hclust(dd)
cor_t_layer_toPlot = cor_t_layer[hc$order, c(1, 7:2)]
rownames(cor_t_layer_toPlot) = ct$lab[match(rownames(cor_t_layer_toPlot), ct$prelimCluster)]
colnames(cor_t_layer_toPlot) = gsub("ayer", "", colnames(cor_t_layer_toPlot))

pdf("pdf/dlpfc_snRNAseq_overlap_heatmap.pdf", width = 10)
print(
    levelplot(
        cor_t_layer_toPlot,
        aspect = "fill",
        at = theSeq,
        col.regions = my.col,
        ylab = "",
        xlab = "",
        scales = list(x = list(rot = 90, cex = 1.5), y = list(cex = 1.5))
    )
)
dev.off()

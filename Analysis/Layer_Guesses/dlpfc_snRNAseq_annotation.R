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

########################
### broad clusters #####
########################

## get pseudobulk
sce.dlpfc$PseudoSample = paste0(sce.dlpfc$sample, ":", sce.dlpfc$prelimCluster)

cIndexes = splitit(sce.dlpfc$PseudoSample)
umiComb <- sapply(cIndexes, function(ii)
    rowSums(assays(sce.dlpfc)$counts[, ii, drop = FALSE]))

phenoComb = colData(sce.dlpfc)[!duplicated(sce.dlpfc$PseudoSample), 13:15]
rownames(phenoComb) = phenoComb$PseudoSample
phenoComb = phenoComb[colnames(umiComb), ]
phenoComb = DataFrame(phenoComb)

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
# 1    1646         770         557
# 2     346         147         103
# 3     315         153         127
# 4     306         118          66
# 5     846         288         165
# 6     247          92          48
# 7     965         281         170
# 8     308         183         142
# 9     302         136          94
# 10    381         159          96
# 11    218          91          51
# 12    356         174         123
# 13    307         177         129
# 14    292         133         104
# 15    207          96          48
# 16    312         114          65
# 17    251         179          31
# 18    198          94          65
# 19    426         149          98
# 20    192          57          17
# 21    173          56          27
# 22    229          86          39
# 23    192          69          39
# 24    609         177         105
# 25    408         153          84
# 26    236          71          28
# 27    274         111          63
# 28    261         119          81
# 29    178          64          33
# 30    318         106          60
# 31    396         157          96

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

dd = dist(1-cor_t_layer)
hc = hclust(dd)
cor_t_layer_toPlot = cor_t_layer[hc$order, c(1, 7:2)]

ct = colData(sce_pseudobulk)
ct = ct[!duplicated(sce_pseudobulk$prelimCluster),]
ct = ct[order(ct$collapsedCluster, ct$prelimCluster),]
ct$lab = paste0(ct$prelimCluster, " (", ct$collapsedCluster,")")

cor_t_layer_toPlot = cor_t_layer[as.character(ct$prelimCluster), c(1, 7:2)]
rownames(cor_t_layer_toPlot) = ct$lab

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

pdf("pdf/dlpfc_snRNAseq_overlap_pheatmap.pdf", width = 10)
print(
    pheatmap(
        cor_t_layer_toPlot,
        ylab = "",
        xlab = "",
    )
)
dev.off()

#### gene expression
g = c("SNAP25", "CAMK2A", "GAD2", "SOX11",
	"FOXP2", "PDGFRA", "MBP", "PLP1",
	"AQP4", "GFAP", "CD74")
t0_contrasts_cell_markers = t0_contrasts_cell[g,]

cc_cell_layer = cor(t(t0_contrasts_cell_markers), cor_t_layer)
signif(cc_cell_layer,3)

### heatmap
theSeq2 = seq(-5, 5, by = 0.01)
my.col2 <- colorRampPalette(brewer.pal(7, "RdBu"))(length(theSeq2))
t0_contrasts_cell_markers_plot = t(t0_contrasts_cell_markers)
t0_contrasts_cell_markers_plot[t0_contrasts_cell_markers_plot > 5] = 5
t0_contrasts_cell_markers_plot = t0_contrasts_cell_markers_plot[,ncol(t0_contrasts_cell_markers_plot):1]

pdf("pdf/dlpfc_snRNAseq_marker_heatmap.pdf", width = 10)
print(
    levelplot(
        t0_contrasts_cell_markers_plot,
        aspect = "fill",
        at = theSeq2,
        col.regions = my.col2,
        ylab = "",
        xlab = "",
        scales = list(x = list(rot = 90, cex = 1.5), y = list(cex = 1.5))
    )
)
dev.off()

########################################
#### specific /collapsed clusters ######
########################################
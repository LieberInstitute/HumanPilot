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
library('limma')

load("rda/sce_layer.Rdata")
source("layer_specificity_functions.R")

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

##############
# load twas ##
##############

##
load("/dcl01/ajaffe/data/lab/dg_hippo_paper/rdas/tt_objects_gene.Rdata")

## split by region
ttList = split(tt, tt$region)
names(ttList)[1] = "DG-GCL"
ttList = lapply(ttList, function(x) {
    x = as.data.frame(x)
    rownames(x) = ss(x$geneid, "\\.")
    x
})

tt_dlpfc = ttList$DLPFC

## compare
mm = match(rownames(tt_dlpfc), rownames(t0_contrasts))
tt_dlpfc_match = tt_dlpfc[!is.na(mm),]
t0_contrasts_match = as.data.frame(t0_contrasts[mm[!is.na(mm)],])
pvals0_contrasts_match = as.data.frame(pvals0_contrasts[mm[!is.na(mm)],])
fdrs0_contrasts_match = as.data.frame(fdrs0_contrasts[mm[!is.na(mm)],])

cor(tt_dlpfc_match$TWAS.Z, t0_contrasts_match)
cor(tt_dlpfc_match$SCZD_t, t0_contrasts_match)
cor(abs(tt_dlpfc_match$SCZD_t), -log10(pvals0_contrasts_match))

## TWAS
layer_twas_enrich_tab = mapply(function(t, f) {
    layer = t > 0 & f < 0.1
    tt = table(TWAS_sig = tt_dlpfc_match$TWAS.Bonf < 0.1, Layer = layer)
},
    t0_contrasts_match,
    fdrs0_contrasts_match,
    SIMPLIFY = FALSE)

fisherList = lapply(layer_twas_enrich_tab, fisher.test, alternative = "greater")
sapply(fisherList, "[[", "p.value")

### cluster the genes in the TWAS significant list
t0_contrasts_TWAS_FDR = as.matrix(t0_contrasts_match[tt_dlpfc_match$TWAS.FDR < 0.05,])
rownames(t0_contrasts_TWAS_FDR) = rowData(sce_layer)$gene_name[match(rownames(t0_contrasts_TWAS_FDR), rownames(sce_layer))]

t0_contrasts_TWAS_Bonf = as.matrix(t0_contrasts_match[tt_dlpfc_match$TWAS.Bonf < 0.1,])
rownames(t0_contrasts_TWAS_Bonf) = rowData(sce_layer)$gene_name[match(rownames(t0_contrasts_TWAS_Bonf), rownames(sce_layer))]

pdf("pdf/TWAS_Heatmaps_FDR.pdf", h = 64, w = 6)
pheatmap(
    t0_contrasts_TWAS_FDR,
    annotation_colors = ann_colors,
    annotation_names_row = TRUE
)
dev.off()
pdf("pdf/TWAS_Heatmaps_Bonf.pdf", h = 20, w = 6)
pheatmap(
    t0_contrasts_TWAS_Bonf,
    annotation_colors = ann_colors,
    annotation_names_row = TRUE
)
dev.off()

##########################
## DE stats

load(
    "/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_noHGoldQSV_matchDLPFC.rda"
)
outGene2 = outGene
outGene2$gencodeTx = NULL

mm2 = match(outGene2$ensemblID, rownames(t0_contrasts))

outGene2_match2 = outGene2[!is.na(mm2),]
t0_contrasts_match2 = as.data.frame(t0_contrasts[mm2[!is.na(mm2)],])
pvals0_contrasts_match2 = as.data.frame(pvals0_contrasts[mm2[!is.na(mm2)],])
fdrs0_contrasts_match2 = as.data.frame(fdrs0_contrasts[mm2[!is.na(mm2)],])

cc_bs2 = cor(outGene2_match2$t, t0_contrasts_match2)
plot(outGene2_match2$t, t0_contrasts_match2$Layer2)
plot(outGene2_match2$t, t0_contrasts_match2$Layer3)
plot(outGene2_match2$t, t0_contrasts_match2$WM)
plot(outGene2_match2$t, t0_contrasts_match2$Layer4)

plot(density(outGene2_match2$t[fdrs0_contrasts_match2$WM > 0.05]))
lines(density(outGene2_match2$t[fdrs0_contrasts_match2$WM < 0.05]), col =
        "red")
plot(density(outGene2_match2$t[fdrs0_contrasts_match2$Layer2 > 0.05]))
lines(density(outGene2_match2$t[fdrs0_contrasts_match2$Layer2 < 0.05]), col =
        "red")
plot(density(outGene2_match2$t[fdrs0_contrasts_match2$Layer3 > 0.05]))
lines(density(outGene2_match2$t[fdrs0_contrasts_match2$Layer3 < 0.05]), col =
        "red")

## DE
layer_de_enrich_tab = mapply(function(t, f) {
    layer = t > 0 & f < 0.1
    tt = table(DE_sig = outGene2_match2$adj.P.Val < 0.05, Layer = layer)
},
    t0_contrasts_match2,
    fdrs0_contrasts_match2,
    SIMPLIFY = FALSE)

fisherList = lapply(layer_de_enrich_tab, fisher.test, alternative = "greater")
bs2 = data.frame(
    P = sapply(fisherList, "[[", "p.value"),
    OR = sapply(layer_de_enrich_tab, getOR)
)

## make heatmap
t0_contrasts_DE_FDR05 = as.matrix(t0_contrasts_match2[outGene2_match2$adj.P.Val < 0.05,])
t0_contrasts_DE_FDR10 = as.matrix(t0_contrasts_match2[outGene2_match2$adj.P.Val < 0.10,])
rownames(t0_contrasts_DE_FDR05) = rowData(sce_layer)$gene_name[match(rownames(t0_contrasts_DE_FDR05), rownames(sce_layer))]
rownames(t0_contrasts_DE_FDR10) = rowData(sce_layer)$gene_name[match(rownames(t0_contrasts_DE_FDR10), rownames(sce_layer))]


pdf("pdf/DE_Heatmaps_FDR05.pdf", h = 24, w = 6)
pheatmap(
    t0_contrasts_DE_FDR05,
    annotation_colors = ann_colors,
    annotation_names_row = TRUE
)
dev.off()

############
## BS1 #####
############
load(
    "/users/ajaffe/Lieber/Projects/RNAseq/SzControl_DE_paper/rdas/DE_statistics_adjAndQsva.rda",
    verbose = TRUE
)
outGene1 = outGene

mm2a = match(names(outGene1), rownames(t0_contrasts))

outGene1_match2 = outGene1[!is.na(mm2a),]
t0_contrasts_match2a = as.data.frame(t0_contrasts[mm2a[!is.na(mm2a)],])
pvals0_contrasts_match2a = as.data.frame(pvals0_contrasts[mm2a[!is.na(mm2a)],])
fdrs0_contrasts_match2a = as.data.frame(fdrs0_contrasts[mm2a[!is.na(mm2a)],])

outGene1_match2 = as.data.frame(outGene1_match2)

cc_bs1 = cor(outGene1_match2[, grep("^tstat", colnames(outGene1_match2))],
    t0_contrasts_match2a)

## DE
layer_de_enrich_tab_adj = mapply(function(t, f) {
    layer = t > 0 & f < 0.1
    tt = table(DE_sig = outGene1_match2$fdr_adj < 0.05, Layer = layer)
},
    t0_contrasts_match2a,
    fdrs0_contrasts_match2a,
    SIMPLIFY = FALSE)
layer_de_enrich_tab_qsva = mapply(function(t, f) {
    layer = t > 0 & f < 0.1
    tt = table(DE_sig = outGene1_match2$fdr_qsva < 0.05, Layer = layer)
},
    t0_contrasts_match2a,
    fdrs0_contrasts_match2a,
    SIMPLIFY = FALSE)
layer_de_enrich_tab_pca = mapply(function(t, f) {
    layer = t > 0 & f < 0.1
    tt = table(DE_sig = outGene1_match2$fdr_pca < 0.1, Layer = layer)
},
    t0_contrasts_match2a,
    fdrs0_contrasts_match2a,
    SIMPLIFY = FALSE)

fisherList_adj = lapply(layer_de_enrich_tab_adj, fisher.test, alternative =
        "greater")
fisherList_qsva = lapply(layer_de_enrich_tab_qsva, fisher.test, alternative =
        "greater")
fisherList_pca = lapply(layer_de_enrich_tab_pca, fisher.test, alternative =
        "greater")

bs1_adj = data.frame(
    P = sapply(fisherList_adj, "[[", "p.value"),
    OR = sapply(layer_de_enrich_tab_adj, getOR)
)
bs1_qsva = data.frame(
    P = sapply(fisherList_qsva, "[[", "p.value"),
    OR = sapply(layer_de_enrich_tab_qsva, getOR)
)
bs1_pca = data.frame(
    P = sapply(fisherList_pca, "[[", "p.value"),
    OR = sapply(layer_de_enrich_tab_pca, getOR)
)

## make heatmap
t0_contrasts_DE_FDR05 = as.matrix(t0_contrasts_match2a[outGene1_match2$adj.P.Val < 0.05,])
t0_contrasts_DE_FDR10 = as.matrix(t0_contrasts_match2a[outGene1_match2$adj.P.Val < 0.10,])
rownames(t0_contrasts_DE_FDR05) = rowData(sce_layer)$gene_name[match(rownames(t0_contrasts_DE_FDR05), rownames(sce_layer))]
rownames(t0_contrasts_DE_FDR10) = rowData(sce_layer)$gene_name[match(rownames(t0_contrasts_DE_FDR10), rownames(sce_layer))]


pdf("pdf/DE_Heatmaps_FDR05.pdf", h = 24, w = 6)
pheatmap(
    t0_contrasts_DE_FDR05,
    annotation_colors = ann_colors,
    annotation_names_row = TRUE
)
dev.off()


##### hg38 rerun
load(
    "/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase1_qsv/rdas/dxStats_dlpfc_filtered_qSVA_BSP1_DLPFC.rda",
    verbose = TRUE
)
outGene1_hg38 = outGene

mm2b = match(outGene1_hg38$ensemblID, rownames(t0_contrasts))

outGene1_match2b = outGene1_hg38[!is.na(mm2b),]
t0_contrasts_match2b = as.data.frame(t0_contrasts[mm2b[!is.na(mm2b)],])
pvals0_contrasts_match2b = as.data.frame(pvals0_contrasts[mm2b[!is.na(mm2b)],])
fdrs0_contrasts_match2b = as.data.frame(fdrs0_contrasts[mm2b[!is.na(mm2b)],])

outGene1_match2b = as.data.frame(outGene1_match2b)

cc_bs1_hg38 = cor(outGene1_match2b$t, t0_contrasts_match2b)

## DE
layer_de_enrich_tab_hg38 = mapply(function(t, f) {
    layer = t > 0 & f < 0.1
    tt = table(DE_sig = outGene1_match2b$adj.P.Val < 0.05, Layer = layer)
},
    t0_contrasts_match2b,
    fdrs0_contrasts_match2b,
    SIMPLIFY = FALSE)

fisherList_hg38 = lapply(layer_de_enrich_tab_hg38, fisher.test, alternative =
        "greater")
bs1_hg38 = data.frame(
    P = sapply(fisherList_hg38, "[[", "p.value"),
    OR = sapply(layer_de_enrich_tab_hg38, getOR)
)

######################
## from psychENCODE ##

library(readxl)

stats = as.data.frame(read_excel("gene_sets/aat8127_Table_S1.xlsx", sheet = "DGE"))
rownames(stats) = stats$ensembl_gene_id

mm3 = match(stats$ensembl_gene_id, rownames(t0_contrasts))
table(is.na(mm3))

stats_match3 = stats[!is.na(mm3),]
t0_contrasts_match3 = as.data.frame(t0_contrasts[mm3[!is.na(mm3)],])
pvals0_contrasts_match3 = as.data.frame(pvals0_contrasts[mm3[!is.na(mm3)],])
fdrs0_contrasts_match3 = as.data.frame(fdrs0_contrasts[mm3[!is.na(mm3)],])

## hmmm
cc_pe = cor(stats_match3[, grep("t.value", colnames(stats_match3))],
    t0_contrasts_match3, use = "pair")

## DE
layer_sczd_enrich_tab = mapply(function(t, f) {
    layer = t > 0 & f < 0.1
    tt = table(DE_sig = stats_match3$SCZ.fdr < 0.05, Layer = layer)
},
    t0_contrasts_match3,
    fdrs0_contrasts_match3,
    SIMPLIFY = FALSE)
fisherList_sczd = lapply(layer_sczd_enrich_tab, fisher.test, alternative =
        "greater")
pe_sczd = data.frame(
    P = sapply(fisherList_sczd, "[[", "p.value"),
    OR = sapply(layer_sczd_enrich_tab, getOR)
)

layer_asd_enrich_tab = mapply(function(t, f) {
    layer = t > 0 & f < 0.1
    tt = table(DE_sig = stats_match3$ASD.fdr < 0.05, Layer = layer)
},
    t0_contrasts_match3,
    fdrs0_contrasts_match3,
    SIMPLIFY = FALSE)
fisherList_asd = lapply(layer_asd_enrich_tab, fisher.test, alternative =
        "greater")
pe_asd = data.frame(
    P = sapply(fisherList_asd, "[[", "p.value"),
    OR = sapply(layer_asd_enrich_tab, getOR)
)

layer_bpd_enrich_tab = mapply(function(t, f) {
    layer = t > 0 & f < 0.1
    tt = table(DE_sig = stats_match3$BD.fdr < 0.05, Layer = layer)
},
    t0_contrasts_match3,
    fdrs0_contrasts_match3,
    SIMPLIFY = FALSE)
fisherList_bpd = lapply(layer_bpd_enrich_tab, fisher.test, alternative =
        "greater")
pe_bpd = data.frame(
    P = sapply(fisherList_bpd, "[[", "p.value"),
    OR = sapply(layer_bpd_enrich_tab, getOR)
)


###############
## combine all DE stats

de_enrich_list = list(
    BS2_SCZD = bs2,
    BS1_SCZD = bs1_qsva,
    BS1_SCZD_remap = bs1_hg38,
    PE_SCZD = pe_sczd,
    PE_ASD = pe_asd,
    PE_BD = pe_bpd
)
de_enrich_tab = do.call("cbind", de_enrich_list)
de_enrich_tab = cbind(de_enrich_tab, t(rbind(cc_bs1, cc_pe)))
colnames(de_enrich_tab)[13:16] = c("BS2.SCZD_Corr", "PE.SCZD_Corr",
    "PE.ASD_Corr", "PE.BD_Corr")

de_enrich_tab = de_enrich_tab[c(2:7, 1),] # make WM bottom
signif(de_enrich_tab, 3)

## make 3 different heatmaps


## TWAS
twas_sczd = as.data.frame(read_excel("gene_sets/aat8127_Table_S4.xlsx", sheet = "SCZ.TWAS"))
rownames(twas_sczd) = twas_sczd$GeneID
twas_sczd_match = twas_sczd[match(rownames(pvals0_contrasts), rownames(twas_sczd)),]
cor(twas_sczd_match$TWAS.Z, t0_contrasts, use = "pair")

twas_asd = as.data.frame(read_excel("gene_sets/aat8127_Table_S4.xlsx", sheet = "ASD.TWAS"))
rownames(twas_asd) = twas_asd$GeneID
twas_asd_match = twas_asd[match(rownames(pvals0_contrasts), rownames(twas_asd)),]
cor(twas_asd_match$TWAS.Z, t0_contrasts, use = "pair")


###############
## images #####
###############


## panel plot
legIndex = c(1, 12, 23, 34, 45, 56)
pdf(
    "publicData/LIBDdlpfc_Stage_markerList_panelPlot.pdf",
    h = 20,
    w = 13
)
gIndexes = splitit(rse_geneDLPFC$ageGroup)
cols = colorRampPalette(brewer.pal(9, "Blues"))(1000)
for (i in seq(along = gIndexes)) {
    ii = gIndexes[[i]]
    yy = rowMeans(yExprsDLPFC_subset[, ii])
    par(mfrow = c(6, 11))
    for (j in 1:nrow(geneTab)) {
        if (j %in% legIndex)
            par(mar = c(1, 4, 3, 0))
        else
            par(mar = c(1, 2, 3, 2))
        layer = c(
            geneTab[j, "1"],
            rep(geneTab[j, "2"], 3),
            rep(geneTab[j, "3"], 2),
            rep(geneTab[j, "4"], 3),
            rep(geneTab[j, "5a"], 3),
            rep(geneTab[j, "5b"], 3),
            rep(geneTab[j, "6"], 8),
            geneTab[j, "SP"]
        )
        layer = layer * yy[j] # weight by expression
        layer[layer > 10] = 10
        layerMat = matrix(rev(layer), nc = 1)
        image(
            t(layerMat),
            col = cols,
            zlim = c(0, 10),
            axes = FALSE,
            main = geneTab$MouseSym[j],
            cex.main = 1.6
        )
        if (j %in% legIndex)
            axis(
                2,
                at = c(0.5, 5, 10.5, 13.5, 16.5, 20, 23.75) / 24,
                labels = c("SP", "VI", "Vb", "Va", "IV", "II/III", "I"),
                cex.axis = 2,
                las = 2
            )
        if (j == 1)
            text(
                x = 0.5,
                y = 0.5,
                names(gIndexes)[i],
                srt = 90,
                cex = 2.5,
                offset = 0
            )
    }
}
dev.off()

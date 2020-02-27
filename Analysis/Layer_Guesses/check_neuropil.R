
####################
# load spot level ##
####################

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
library('janitor')
library('org.Hs.eg.db')
library('GenomicFeatures')
library('SparseM')
library('biomaRt')
library('rtracklayer')
library('liftOver')
library('readr')
library('R.utils')
library('RColorBrewer')
library('spatialLIBD')
library('M3C')

## Load data
load(here(
    'Analysis',
    'Human_DLPFC_Visium_processedData_sce_scran.Rdata'
))

## For plotting
source(here('Analysis', 'spatialLIBD_global_plot_code.R'))
genes <- paste0(rowData(sce)$gene_name, '; ', rowData(sce)$gene_id)

## Functions derived from this script, to make it easier to resume the work
sce_layer_file <-
    here('Analysis', 'Layer_Guesses', 'rda', 'sce_layer.Rdata')
if (file.exists(sce_layer_file))
    load(sce_layer_file, verbose = TRUE)
source(here('Analysis', 'Layer_Guesses', 'layer_specificity_functions.R'))

## Load layer guesses
load(here('Analysis', 'Layer_Guesses','rda',
    'layer_guess_tab.Rdata'))

## Add layer guesses to the sce object
sce$layer_guess <- NA
m <- match(sce$key, layer_guess_tab$key)
table(is.na(m))
# FALSE  TRUE
# 47329   352
sce$layer_guess[!is.na(m)] <- layer_guess_tab$layer[m[!is.na(m)]]

## Drop the layer guess NAs for now
sce_original <- sce
sce <- sce[, !is.na(sce$layer_guess)]
dim(sce)
# [1] 33538 47329

## Next, re-label "Layer 2/3" as "Layer 3" for now
## (there's more layer 3 in the other samples than 2 anyway)
sce$layer_guess[sce$layer_guess == 'Layer 2/3'] <- 'Layer 3'

## Make it into a factor with WM as the reference
## and remove spaces
sce$layer_guess <-
    factor(gsub(' ', '', sce$layer_guess), levels = c('WM', paste0('Layer', 1:6)))

#############################
### line up to layer info ###
#############################

## get out counts

## filter
exprsIndex = rowMeans(assays(sce)$logcounts) > 0
sce = sce[exprsIndex,]

## fit model of 0 cell bodies vs all others
sce$zero_cell = as.numeric(sce$cell_count  == 0)

###############################
##### EDA on Neuropil spots ###
###############################
sce0 = sce[,sce$zero_cell == 1]
sce0_hvg_logcounts = as.matrix(assays(sce0)$logcounts[top.hvgs,])
sce0_logcounts = as.matrix(assays(sce0)$logcounts)

g_pre = toupper(c("Bsn", "Rims1", "Rims2", "Rims3", "Stx6", "Rapgef4"))
g_post = toupper(c("Arc", "Bdnf", "Grin1", "SLC17A7"))
g_other = toupper(c("Dlg4", "Grin2a", "Limk1"))
g_ecm = toupper(c("HAPLN1", "ACAN", "BCAN", "TNR"))
g_my = toupper(c("Mobp", "Mbp", "Aqp4"))
g = c(g_pre, g_post,g_other,g_ecm,g_my)
g %in% rowData(sce0)$gene_name

sce0_logcounts_syn = sce0_logcounts[match(g, rowData(sce0)$gene_name),]
rownames(sce0_logcounts_syn) = g

pca_marker = prcomp(t(sce0_logcounts_syn))
pcaVars_marker = getPcaVars(pca_marker)

palette(libd_layer_colors[1:7])
plot(pca_marker$x, pch = 21, bg = sce0$layer_guess)
round(pca_marker$rot[,1:3],3)

plot(pca_marker$x[,3:2], pch = 21, bg = sce0$layer_guess)

plot(pca_marker$x[,2:3], pch = 21, bg = sce0$layer_guess)
plot(pca_marker$x[,c(1,3)], pch = 21, bg = sce0$layer_guess)

### 
tsne(t(sce0_logcounts_syn))

dd = dist(t(sce0_logcounts_syn))
hc = hclust(dd)
myplclust(hc, lab.col = as.numeric(sce0$layer_guess),cex=0.33)
table(cutree(hc, h=8.25))

g = cutree(hc, h=8.25) 
g[g==7] = NA
palette(brewer.pal(6,"Dark2"))
myplclust(hc, lab.col = g,cex=0.33)

plot(pca0$x[,1], sce0$sum_umi, pch = 21, bg = sce0$layer_guess)
legend("bottomright", levels(sce0$layer_guess), col = 1:7, pch = 15)
cc = cor(sce0_hvg_logcounts)


#########################
### read in neuropil ####
#########################

######## old paper
# neuropil = read_excel("gene_sets/1-s2.0-S0896627312002863-mmc2.xlsx", sheet=1, skip=29)
# neuropil = clean_names(neuropil)
# neuropil = as.data.frame(neuropil)

## homology
# # # hom = read.delim("http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt",header=FALSE,
# hom = read.delim("http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt",
	# as.is=TRUE)


# hom_hs = hom[hom$Common.Organism.Name == "human",]
# neuropil$entrez_hs = hom_hs$EntrezGene.ID[match(neuropil$HomoloGene.ID, hom_hs$HomoloGene.ID)]
# neuropil$symbol_hs = hom_hs$Symbol[match(neuropil$HomoloGene.ID, hom_hs$HomoloGene.ID)]

# hom_rat = hom[match(neuropil$entrez, hom$EntrezGene.ID),]
# hom_hs = hom[hom$Common.Organism.Name == "human",]
# neuropil$entrez_hs = hom_hs$EntrezGene.ID[match(hom_rat$HomoloGene.ID, hom_hs$HomoloGene.ID)]
# neuropil$symbol_hs = hom_hs$Symbol[match(hom_rat$HomoloGene.ID, hom_hs$HomoloGene.ID)]

##### new paper
neuropil = read_excel("gene_sets/aau3644_TableS1.xls", 
	sheet="DE", skip=4)
neuropil = as.data.frame(neuropil)
neuropil$padj = as.numeric(neuropil$padj)
table(neuropil$padj < 0.05, sign(neuropil$stat))

## line up
hom2 = read.delim("http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt",
	as.is=TRUE, header=FALSE)
colnames(hom2)[1:5] = c("symbol_hs", "entrez_hs", "entrez_mm", "mapped", "symbol_mm")
neuropil$symbol_hs = hom2$symbol_hs[match(neuropil$GeneID, hom2$symbol_mm)]
neuropil$entrez_hs = hom2$entrez_hs[match(neuropil$GeneID, hom2$symbol_mm)]

# get ensembl from entrez
ens = select(org.Hs.eg.db, key=as.character(neuropil$entrez_hs),
	columns="ENSEMBL", keytype="ENTREZID")
neuropil$ensembl_hs = ens$ENSEMBL[match(neuropil$entrez_hs, ens$ENTREZID)]

# #### my stats recaluclated
# neuropil_2 = read.csv("../hafner_vglut/vglut1_de_stats.csv",
	# row.names=1, as.is=TRUE)
# neuropil_2$symbol_hs = hom2$symbol_hs[match(neuropil_2$EntrezID, hom2$entrez_mm)]
# neuropil_2$entrez_hs = hom2$entrez_hs[match(neuropil_2$EntrezID, hom2$entrez_mm)]
# neuropil_2$ensembl_hs = ens$ENSEMBL[match(neuropil_2$entrez_hs, ens$ENTREZID)]

# ensembl=useMart("ensembl",dataset = 'mmusculus_gene_ensembl')
# MMtoHG = getBM(attributes = c('ensembl_gene_id','hsapiens_homolog_ensembl_gene'),mart = ensembl)
# neuropil_2$ensembl_hs = MMtoHG$hsapiens_homolog_ensembl_gene[match(

#######################################
### spot level below, very slow, run once and save
sce_logcounts = assays(sce)$logcounts

# mod = model.matrix(~zero_cell + layer_guess + subject_position,data =colData(sce))
# fit = bumphunter:::.getEstimate(sce_logcounts, mod, coef=2, full=TRUE)
# tt = bumphunter:::.getModT(fit)
# pv = 2*pt(-abs(tt$t), df = tt$df.total)
# stats = data.frame(logFC = fit$coef, t = tt$t, P.Value = pv)
# rownames(stats) = rownames(sce_logcounts)
# stats$Symbol = rowData(sce)$gene_name

# ## by layer
# layerIdx = splitit(sce$layer_guess)
# statList_layer = mclapply(layerIdx,  function(ii) {
	# mod = model.matrix(~zero_cell + subject_position,
		# data =colData(sce[,ii]))
	# fit = bumphunter:::.getEstimate(
		# sce_logcounts[,ii], mod, coef=2, full=TRUE)
	# tt = bumphunter:::.getModT(fit)
	# pv = 2*pt(-abs(tt$t), df = tt$df.total)
	# stats = data.frame(logFC = fit$coef, t = tt$t, P.Value = pv)
	# rownames(stats) = rownames(sce_logcounts)
	# stats
# }, mc.cores=7)

# save(stats, statList_layer,
	# file = "rda/linear_model_zerospot_byLayer.Rdata")

## human anno
load("rda/linear_model_zerospot_byLayer.Rdata")

write.csv(stats, file = "SupplementaryTableXX_spot_stats.csv")
stats$in_neuropil = rownames(stats) %in% neuropil$ensembl_hs[neuropil$padj < 0.05]
stats$neuropil_padj = neuropil$padj[match(rownames(stats), neuropil$ensembl_hs)]
stats$neuropil_stat = neuropil$stat[match(rownames(stats), neuropil$ensembl_hs)]
stats$neuropil_logFC = neuropil$Log2FoldChange[match(rownames(stats), neuropil$ensembl_hs)]

## plots
ct = cor.test(stats$logFC, stats$neuropil_logFC)
ct$p.value # 8.002645e-160
ct_sig = with(stats[which(stats$neuropil_padj < 0.05),], cor.test(logFC, neuropil_logFC))
ct_sig$p.value # 1.898373e-30

n_all = sum(complete.cases(stats[,c("logFC", "neuropil_logFC")]))
n_sig = sum(complete.cases(stats[which(stats$neuropil_padj < 0.05),c("logFC", "neuropil_logFC")]))

pdf("pdf/neuropil_scatter_plots.pdf")
par(mar=c(5,6,4,2),cex.axis=1.8,cex.lab=1.8,cex.main=1.8)
plot(logFC ~ neuropil_logFC, data=stats,
	pch = 21, bg = "grey", 
	main = paste(n_all, "Expressed Homologs"),
	xlab = "vGLUT1+ Enrichment (log2FC)",
	ylab = "Visium Zero Cell Enrichment (log2FC)")
legend("topleft", paste0("r=", signif(ct$estimate, 3)),cex=1.5)
plot(logFC ~ neuropil_logFC, 
	data=stats[which(stats$neuropil_padj < 0.05),],
	main = paste(n_sig, "vGLUT1+ Significant Homologs"),
	pch = 21, bg = "grey",
	xlab = "vGLUT1+ Enrichment (log2FC)",
	ylab = "Visium Zero Cell Enrichment (logFC)")
legend("topleft", paste0("r=", signif(ct_sig$estimate, 3)),cex=1.5)
dev.off()

### check our TWAS
load("/dcl01/ajaffe/data/lab/dg_hippo_paper/rdas/tt_objects_gene.Rdata")
tt_dlpfc=  as.data.frame(tt[tt$region == "DLPFC",])
tt_dlpfc$ensemblID = ss(tt_dlpfc$geneid, "\\.")

tt_dlpfc_2 = tt_dlpfc[match(rownames(stats), ss(tt_dlpfc$geneid, "\\.")),]
		
plot(-log10(stats$P.Value), -log10(tt_dlpfc_2$TWAS.P),
	xlim = c(0,20))
	
cor.test(stats$t, abs(tt_dlpfc_2$TWAS.Z))
lapply(statList_layer, function(x) cor.test(x$t, abs(tt_dlpfc_2$TWAS.Z)))

##################
### make LDSC gene lists
#################

statList_layer$All = stats
statList = lapply(statList_layer, function(x) {
	x$adj.P.Val = p.adjust(x$P.Value,"fdr")
	x 
})

## get out genes
geneList_enr = lapply(statList, function(x) {
	rownames(x)[x$adj.P.Val < 0.1 & x$t > 0]
})
geneList_dep = lapply(statList, function(x) {
	rownames(x)[x$adj.P.Val < 0.1 & x$t < 0]
})

## combine
names(geneList_enr) = paste0(names(geneList_enr), "_Enriched")
names(geneList_dep) = paste0(names(geneList_dep), "_Depleted")
geneList = c(geneList_enr, geneList_dep)
lengths(geneList)


## get coordinates, hg19 for LDSC
gtf = import("/dcl01/ajaffe/data/lab/singleCell/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf")
gtf = gtf[gtf$type == "gene"]
seqlevels(gtf)[1:25] = paste0("chr", seqlevels(gtf)[1:25])
names(gtf) = gtf$gene_id
gtf = gtf[rownames(stats)]

## we should add some random gene lists
set.seed(34135)
n = rep(c(10,50,100,200,300,500,1000,2000,3000), each=3)
grListNull = GRangesList(lapply(n, function(nn) 
	granges(gtf[sample(1:length(gtf), nn)])))
names(grListNull) = paste0("null_", n, "_rep", rep(1:3, times = 9))

## list of GRanges for all real genes
grList = GRangesList(lapply(geneList, function(g) granges(gtf[g])))

## combine
grList = c(grList, grListNull)

## pad with promoter sequence
grList = endoapply(grList, promoters, upstream = 5000,
		downstream = 2000,  use.names=TRUE)

## lift to hg19
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)
seqlevelsStyle(grList) = "UCSC"  # necessary

categories = endoapply(grList,function(g) {
	gg = unlist(liftOver(g,ch))
	genome(gg) = "hg19"
	gg
})

saveRDS(categories, "LDSC/categories.rds")
write.table(data.frame(names(categories)), 
	file = "LDSC/categories.txt",
	quote=FALSE, row.names=FALSE,col.names=FALSE)


### ============================================================================
### Make annotations
###
seqlevels <- 1:22

# NOTE: Don't re-make the CNS annotation
k <- which(names(categories)!="CNS")

mclapply(seqlevels, function(sl) {
  message(sl)
  cds <- read_tsv(
    paste0("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1/cell_type_groups/CNS.", sl,
           ".annot.gz"))
  cds_gr <- GRanges(paste0("chr", cds$CHR), IRanges(cds$BP, width = 1L))
  annot <- cds[, c("CHR", "BP", "SNP", "CM")]
  
  ol = sapply(categories, countOverlaps, query = cds_gr)
  ol[ol > 1] = 1
  annot = cbind(annot, ol)
  
  # 'Marginal' annotation file
  mclapply(names(categories), function(n) {
    fl <- paste0("LDSC/LDScore/", n, ".Phase1.", sl, ".annot")
    write_tsv(annot[, c("CHR", "BP", "SNP", "CM", n)], fl)
    gzip(fl, overwrite = TRUE)
  }, mc.cores = 4)
}, mc.cores = 10)


# ------------------------------------------------------------------------------
# Annot file for 'base' category
# NOTE: Not added to 'categories' since it is not a brain-derived category
#

mclapply(seqlevels, function(sl) {
  x <- read_tsv(file.path("/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/", "Phase1", "baseline",
                          paste0("baseline.", sl, ".annot.gz")))
  write_tsv(x[, 1:5], paste0("/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/LDSC/LDScore/base.Phase1.", sl, ".annot.gz"))
}, mc.cores = 4)


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
library('MatrixModels')
library('biomaRt')

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

#### my stats recaluclated
neuropil_2 = read.csv("../hafner_vglut/vglut1_de_stats.csv",
	row.names=1, as.is=TRUE)
neuropil_2$symbol_hs = hom2$symbol_hs[match(neuropil_2$EntrezID, hom2$entrez_mm)]
neuropil_2$entrez_hs = hom2$entrez_hs[match(neuropil_2$EntrezID, hom2$entrez_mm)]
neuropil_2$ensembl_hs = ens$ENSEMBL[match(neuropil_2$entrez_hs, ens$ENTREZID)]

# ensembl=useMart("ensembl",dataset = 'mmusculus_gene_ensembl')
# MMtoHG = getBM(attributes = c('ensembl_gene_id','hsapiens_homolog_ensembl_gene'),mart = ensembl)
# neuropil_2$ensembl_hs = MMtoHG$hsapiens_homolog_ensembl_gene[match(

#######################################
### spot level below, very slow, run once and save
# sce_logcounts = assays(sce)$logcounts
# mod = model.matrix(~zero_cell + layer_guess + subject_position,data =colData(sce))
# fit = bumphunter:::.getEstimate(sce_logcounts, mod, coef=2, full=TRUE)
# tt = bumphunter:::.getModT(fit)
# pv = 2*pt(-abs(tt$t), df = tt$df.total)
# stats = data.frame(logFC = fit$coef, t = tt$t, P.Value = pv)
# rownames(stats) = rownames(sce_logcounts)
# stats$Symbol = rowData(sce)$gene_name
# save(stats, file = "rda/linear_model_zerospot.Rdata")

## human anno
load("rda/linear_model_zerospot.Rdata")
stats_2 = stats # for second dataset

stats$in_neuropil = rownames(stats) %in% neuropil$ensembl_hs[neuropil$padj < 0.05]
stats$neuropil_padj = neuropil$padj[match(rownames(stats), neuropil$ensembl_hs)]
stats$neuropil_stat = neuropil$stat[match(rownames(stats), neuropil$ensembl_hs)]
stats$neuropil_logFC = neuropil$Log2FoldChange[match(rownames(stats), neuropil$ensembl_hs)]

stats_2$in_neuropil = rownames(stats) %in% neuropil_2$ensembl_hs[neuropil_2$adj.P.Val < 0.05]
stats_2$neuropil_padj = neuropil_2$adj.P.Val[match(rownames(stats_2), neuropil_2$ensembl_hs)]
stats_2$neuropil_stat = neuropil_2$t[match(rownames(stats_2), neuropil_2$ensembl_hs)]
stats_2$neuropil_logFC = neuropil_2$logFC[match(rownames(stats_2), neuropil_2$ensembl_hs)]

## plots
ct = cor.test(stats$logFC, stats$neuropil_logFC)
ct_sig = with(stats[which(stats$neuropil_padj < 0.05),], cor.test(logFC, neuropil_logFC))
ct_sig$p.value
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
	main = paste(n_sig, "vGLUT1+ Significat Homologs"),
	pch = 21, bg = "grey",
	xlab = "vGLUT1+ Enrichment (log2FC)",
	ylab = "Visium Zero Cell Enrichment (logFC)")
legend("topleft", paste0("r=", signif(ct_sig$estimate, 3)),cex=1.5)
dev.off()

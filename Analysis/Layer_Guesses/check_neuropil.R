
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

stats$in_neuropil = rownames(stats) %in% neuropil$ensembl_hs[neuropil$padj < 0.05]
stats$neuropil_padj = neuropil$padj[match(rownames(stats), neuropil$ensembl_hs)]
stats$neuropil_stat = neuropil$stat[match(rownames(stats), neuropil$ensembl_hs)]
stats$neuropil_logFC = neuropil$Log2FoldChange[match(rownames(stats), neuropil$ensembl_hs)]

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


plot(t ~ neuropil_stat, data=stats)


boxplot(logFC ~ in_neuropil, data=stats)



########################################
## pseudobulk by layer and cell count 
########################################

layerIndexes <-  splitit(paste(sce$sample_name, sce$layer_guess, sce$zero_cell,sep='_'))

## Collapse UMIs
umiComb <-
    sapply(layerIndexes, function(ii)
        rowSums(assays(sce)$counts[, ii, drop = FALSE]))
dim(umiComb)
# [1] 25595   145

## Build a data.frame with the pheno data by layer
layer_df <- data.frame(
    sample_name = ss(colnames(umiComb), '_', 1),
    layer_guess = factor(ss(colnames(umiComb), '_', 2), levels = levels(sce$layer_guess)),
    zero_cells = ss(colnames(umiComb), '_', 3),
    stringsAsFactors = FALSE
)
m_layer <- match(layer_df$sample_name, sce$sample_name)
layer_df$subject <- sce$subject[m_layer]
layer_df$position <- sce$position[m_layer]
layer_df$replicate <- sce$replicate[m_layer]
layer_df$subject_position <- sce$subject_position[m_layer]
rownames(layer_df) <- colnames(umiComb)

## Build a new sce object
sce_layer_count <-
    logNormCounts(SingleCellExperiment(
        list(counts = umiComb),
        colData = layer_df,
        rowData = rowData(sce)
    ))
	

##################
### modeling #####
##################

## Extract the data
mat <- assays(sce_layer_count)$logcounts

## Build a group model
mod <- with(colData(sce_layer_count), model.matrix(~ zero_cells + layer_guess))
colnames(mod) <- gsub('layer_guess', '', colnames(mod))

## Takes like 2 min to run
corfit <-
    duplicateCorrelation(mat, mod, block = sce_layer_count$subject_position)
fit <-
    lmFit(
        mat,
        design = mod,
        block = sce_layer_count$subject_position,
        correlation = corfit$consensus.correlation
    )
eb <- eBayes(fit)

outGene_zeroCell = topTable(eb, coef=2,sort="none", p=1,n=nrow(mat))
outGene_zeroCell$Symbol = rowData(sce_layer_count)$gene_name


hist(outGene_zeroCell$t, xlab="0 Cell Effect");abline(v=0,col="red")
outGene_zeroCell$Symbol[outGene_zeroCell$t > 3]

######################
#### try with voom ###
######################
library(edgeR)
dge = DGEList(counts = assays(sce_layer_count)$counts, 
	genes = rowData(sce_layer_count))
dge = calcNormFactors(dge)

vGene = voom(dge,mod,plot=TRUE)
corfit2 <- duplicateCorrelation(vGene$E, mod, block = sce_layer_count$subject_position)

fitGene = lmFit(vGene, 
	correlation=corfit2$consensus.correlation, 
	block=colData(sce_layer_count)$subject_position)
	
ebGene = eBayes(fitGene)
outGene_zeroCell2 = topTable(ebGene, coef=2, 
	p.value = 1,sort="none",number=nrow(dge))

hist(outGene_zeroCell2$t, xlab="0 Cell Effect");abline(v=0,col="red")
	
###########################
### line up to neuropil data

mm = match(rownames(outGene_zeroCell), neuropil$ensembl_hs)
outGene_zeroCell$in_neuropil = rownames(outGene_zeroCell) %in% neuropil$ensembl_hs[neuropil$padj < 0.05]
outGene_zeroCell$neuropil_stat = neuropil$stat[mm]
outGene_zeroCell$neuropil_logFC = neuropil$Log2FoldChange[mm]
outGene_zeroCell2$in_neuropil = rownames(outGene_zeroCell2) %in% neuropil$ensembl_hs[neuropil$padj < 0.05]
outGene_zeroCell2$neuropil_stat = neuropil$stat[mm]
outGene_zeroCell2$neuropil_logFC = neuropil$Log2FoldChange[mm]

boxplot(t ~ in_neuropil, data=outGene_zeroCell2)
boxplot(t ~ in_neuropil, data=outGene_zeroCell)

plot(t ~ neuropil_stat, data=outGene_zeroCell)
plot(logFC ~ neuropil_logFC, data=outGene_zeroCell)
cor.test(outGene_zeroCell$logFC, outGene_zeroCell$neuropil_logFC)

plot(t ~ neuropil_stat, data=outGene_zeroCell2)
plot(logFC ~ neuropil_logFC, data=outGene_zeroCell2)
cor.test(outGene_zeroCell2$logFC, outGene_zeroCell2$neuropil_logFC)



stats[order(stats$t, decreasing=TRUE)[1:20],]

## this is going to be way too slow, so lets just do linear
# corfit <-  duplicateCorrelation(sce_logcounts, mod, block = sce$subject_position)


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

#########################
### read in neuropil ####
#########################

neuropil = read_excel("gene_sets/1-s2.0-S0896627312002863-mmc2.xlsx", sheet=1, skip=29)
neuropil = clean_names(neuropil)
neuropil = as.data.frame(neuropil)

### homology
# hom = read.delim("http://www.informatics.jax.org/downloads/reports/HMD_HumanPhenotype.rpt",header=FALSE,
hom = read.delim("http://www.informatics.jax.org/downloads/reports/HOM_AllOrganism.rpt",
	as.is=TRUE)
hom_rat = hom[match(neuropil$entrez, hom$EntrezGene.ID),]
hom_hs = hom[hom$Common.Organism.Name == "human",]
neuropil$entrez_hs = hom_hs$EntrezGene.ID[match(hom_rat$HomoloGene.ID, hom_hs$HomoloGene.ID)]
neuropil$symbol_hs = hom_hs$Symbol[match(hom_rat$HomoloGene.ID, hom_hs$HomoloGene.ID)]

## get ensembl from entrez
ens = select(org.Hs.eg.db, key=as.character(neuropil$entrez_hs),
	columns="ENSEMBL", keytype="ENTREZID")
neuropil$ensembl_hs = ens$ENSEMBL[match(neuropil$entrez_hs, ens$ENTREZID)]


#############################
### line up to layer info ###
#############################

## get out counts

## filter
exprsIndex = rowMeans(assays(sce)$logcounts) > 0
sce = sce[exprsIndex,]

## fit model of 0 cell bodies vs all others
sce$zero_cell = as.numeric(sce$cell_count  == 0)

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

## human anno
mm = match(rownames(outGene_zeroCell2), neuropil$ensembl_hs)
neuropil_match = neuropil[mm,]

outGene_zeroCell2$neuropil_numReads = neuropil_match$number_reads_in_cds
outGene_zeroCell2$neuropil_numReads[is.na(outGene_zeroCell2$neuropil_numReads)] = 0

outGene_zeroCell2$neuropil_analyzed = neuropil_match$status == "analyzed"
outGene_zeroCell2$neuropil_analyzed[is.na(outGene_zeroCell2$neuropil_analyzed)] = 0

plot(logFC ~ log2(neuropil_numReads+1), data=outGene_zeroCell2)
cor.test(outGene_zeroCell2$logFC, log2(outGene_zeroCell2$neuropil_numReads+1))

boxplot(t ~ neuropil_analyzed, data=outGene_zeroCell2)
boxplot(t ~ neuropil_analyzed, data=outGene_zeroCell)

#######################################
### spot level below, very slow

mod = model.matrix(~zero_cell + layer_guess + subject_position,data =colData(sce))
fit = bumphunter:::.getEstimate(sce_logcounts, mod, coef=2, full=TRUE)
tt = bumphunter:::.getModT(fit)
pv = 2*pt(-abs(tt$t), df = tt$df.total)
stats = data.frame(logFC = fit$coef, t = tt$t, P.Value = pv)
rownames(stats) = rownames(sce_logcounts)
stats$Symbol = rowData(sce)$gene_name
save(stats, file = "rda/linear_model_zerospot.Rdata")

## human anno
mm2 = match(rownames(stats), neuropil$ensembl_hs)
neuropil_match2 = neuropil[mm2,]

stats$neuropil_numReads = neuropil_match2$number_reads_in_cds
stats$neuropil_numReads[is.na(stats$neuropil_numReads)] = 0

stats$neuropil_analyzed = neuropil_match2$status == "analyzed"
stats$neuropil_analyzed[is.na(stats$neuropil_analyzed)] = 0

plot(logFC ~ log2(neuropil_numReads+1), data=stats)
cor.test(stats$logFC, log2(stats$neuropil_numReads+1))

boxplot(t ~ neuropil_analyzed, data=stats)

stats[order(stats$t, decreasing=TRUE)[1:20],]

## this is going to be way too slow, so lets just do linear
# corfit <-  duplicateCorrelation(sce_logcounts, mod, block = sce$subject_position)

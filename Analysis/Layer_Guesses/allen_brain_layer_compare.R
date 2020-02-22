###
#  module load conda_R/3.6.x
library(jaffelab)
library(Seurat)
library(scater)
library(DropletUtils)
library(limma)
library(scrattch.io)
options(stringsAsFactors = FALSE)
library(org.Hs.eg.db)
library(GenomicFeatures)
library(vroom)
library(Matrix)
library(lattice)
library(RColorBrewer)

## location of data on cluster
## data via: https://portal.brain-map.org/atlases-and-data/rnaseq
## and  https://celltypes.brain-map.org/api/v2/well_known/694416044

path = "/dcl01/ajaffe/data/lab/singleCell/allen_human/data/"

## annotation via: //celltypes.brain-map.org/api/v2/well_known/694416044
gene_map = read.csv(paste0(path, "human_MTG_2018-06-14_genes-rows.csv"), as.is =
        TRUE)
## add ensembl
ens = select(org.Hs.eg.db,
    columns = c("ENSEMBL", "ENTREZID"),
    keys = as.character(unique(gene_map$entrez_id)))
gene_map$ensemblID = ens$ENSEMBL[match(gene_map$entrez_id, ens$ENTREZID)]

##############################
##### smaller MTG 16000 ######
##############################

pheno_mtg = read.csv(
    paste0(path, "human_MTG_2018-06-14_samples-columns.csv"),
    as.is = TRUE,
    row.names = 1
)

## indicate postmortem vs neurosurgery
pheno_mtg$postmortem = ifelse(pheno_mtg$donor %in% c("H16.03.004", "H16.06.002", "H16.06.008", "H16.06.009"),
    0,
    1)

## cell types
pheno_mtg$CellType = ifelse(pheno_mtg$facs_sort_criteria == "NeuN-negative",
    "Glial",
    "Neuronal")

## pseudobulking
pheno_mtg$PseudoSample = with(pheno_mtg, paste(donor, brain_subregion, CellType, sep =
        ":"))

## read in genomics data
exons_mtg = vroom(paste0(path, "human_MTG_2018-06-14_exon-matrix.csv"), delim = ",")
colnames(exons_mtg)[1] = "Gene_ID"
introns_mtg = vroom(paste0(path, "human_MTG_2018-06-14_intron-matrix.csv"),
    delim = ",")
colnames(introns_mtg)[1] = "Gene_ID"

exons_mtg = as.data.frame(exons_mtg)
rownames(exons_mtg) = exons_mtg$Gene_ID
exons_mtg$Gene_ID = NULL
introns_mtg = as.data.frame(introns_mtg)
rownames(introns_mtg) = introns_mtg$Gene_ID
introns_mtg$Gene_ID = NULL

## make sparse
introns_mtg = Matrix(as.matrix(introns_mtg), sparse = TRUE)
exons_mtg = Matrix(as.matrix(exons_mtg), sparse = TRUE)

## combine
total_mtg = introns_mtg + exons_mtg
identical(colnames(total_mtg), rownames(pheno_mtg))
## convert
sce <-   SingleCellExperiment(list(counts = total_mtg),
    colData = pheno_mtg)

rowData(sce) = gene_map[match(rownames(sce), gene_map$entrez_id), ]

save(sce, file = "rda/allen_snRNAseq_sce_MTG.Rdata")

#############
## split and combine
cIndexes = splitit(sce$PseudoSample)
umiComb <- sapply(cIndexes, function(ii)
    rowSums(assays(sce)$counts[, ii, drop = FALSE]))

phenoComb = colData(sce)[!duplicated(sce$PseudoSample), c(4:9, 34:36)]

rownames(phenoComb) = phenoComb$PseudoSample
phenoComb = phenoComb[colnames(umiComb), ]
phenoComb = DataFrame(phenoComb)

sce_pseudobulk <-
    logNormCounts(SingleCellExperiment(
        list(counts = umiComb),
        colData = phenoComb,
        rowData = rowData(sce)
    ))

save(sce_pseudobulk, file = "rda/allen_snRNAseq_sce_MTG_pseudobulked.Rdata")

###############
## analysis ###
###############

load("rda/allen_snRNAseq_sce_MTG_pseudobulked.Rdata")

## drop surgical
sce_pseudobulk = sce_pseudobulk[, sce_pseudobulk$postmortem == 1]
mat <- assays(sce_pseudobulk)$logcounts

## model
mod <- with(colData(sce_pseudobulk),
    model.matrix(~ 0 + brain_subregion + CellType))
colnames(mod) <- gsub('brain_subregion', '', colnames(mod))
colnames(mod) <- gsub('CellType', '', colnames(mod))

## get duplicate correlation
corfit <-
    duplicateCorrelation(mat, mod, block = sce_pseudobulk$donor)
save(corfit, file = "rda/allen_snRNAseq_MTG_pseudobulked_dupCor.Rdata")

## Next for each layer test that layer vs the rest
layer_idx <- splitit(sce_pseudobulk$brain_subregion)

eb0_list_layer <- lapply(layer_idx, function(x) {
    res <- rep(0, ncol(sce_pseudobulk))
    res[x] <- 1
    m <- with(colData(sce_pseudobulk),
        model.matrix(~ res + CellType))
    
    eBayes(
        lmFit(
            mat,
            design = m,
            block = sce_pseudobulk$external_donor_name_label,
            correlation = corfit$consensus.correlation
        )
    )
})
save(eb0_list_layer, file = "rda/allen_snRNAseq_MTG_pseudobulked_specific_Ts.Rdata")

##########
## Extract the p-values
load("rda/allen_snRNAseq_MTG_pseudobulked_specific_Ts.Rdata")

pvals0_contrasts_layer <- sapply(eb0_list_layer, function(x) {
    x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts_layer) = rownames(mat)

t0_contrasts_layer <- sapply(eb0_list_layer, function(x) {
    x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts_layer) = rownames(mat)
fdrs0_contrasts_layer = apply(pvals0_contrasts_layer, 2, p.adjust, 'fdr')

data.frame(
    'FDRsig' = colSums(fdrs0_contrasts_layer < 0.05 &
            t0_contrasts_layer > 0),
    'Pval10-6sig' = colSums(pvals0_contrasts_layer < 1e-6 &
            t0_contrasts_layer > 0),
    'Pval10-8sig' = colSums(pvals0_contrasts_layer < 1e-8 &
            t0_contrasts_layer > 0)
)
# FDRsig Pval10.6sig Pval10.8sig
# L1      4           2           0
# L2      0           0           0
# L3      0           0           0
# L4      0           0           0
# L5      1           1           0
# L6      5           3           0


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

mm = match(rownames(pvals0_contrasts),
    rowData(sce_pseudobulk)$ensemblID)

pvals0_contrasts = pvals0_contrasts[!is.na(mm), ]
t0_contrasts = t0_contrasts[!is.na(mm), ]
fdrs0_contrasts = fdrs0_contrasts[!is.na(mm), ]

pvals0_contrasts_layer = pvals0_contrasts_layer[mm[!is.na(mm)], ]
t0_contrasts_layer = t0_contrasts_layer[mm[!is.na(mm)], ]
fdrs0_contrasts_layer = fdrs0_contrasts_layer[mm[!is.na(mm)], ]

cor_t = cor(t0_contrasts_layer, t0_contrasts)
signif(cor_t, 2)

### just layer specific genes from ones left
layer_specific_indices = mapply(function(t, p) {
    oo = order(t, decreasing = TRUE)[1:100]
},
    as.data.frame(t0_contrasts),
    as.data.frame(pvals0_contrasts))
layer_ind = unique(as.numeric(layer_specific_indices))

cor_t_layer = cor(t0_contrasts_layer[layer_ind, ],
    t0_contrasts[layer_ind, ])
signif(cor_t_layer, 2)


### heatmap
theSeq = seq(-.85, .85, by = 0.01)
my.col <- colorRampPalette(brewer.pal(7, "PRGn"))(length(theSeq))

cor_t_layer_toPlot = cor_t_layer[,c(1, 7:2)]
pdf("pdf/allen_brain_layer_overlap_heatmap.pdf", width = 8)
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


##############################
#### larger 50000 sample #####
##############################

## read in genomic data
tome = paste0(path, "transcrip.tome")
exons = read_tome_dgCMatrix(tome, "data/t_exon")    # (or data/exon)
introns = read_tome_dgCMatrix(tome, "data/t_intron")  # (or data/intron)
total = exons + introns # combine counts
## get annotation info
sample_name = read_tome_sample_names(tome)
colnames(total) = sample_name
gene_name = read_tome_gene_names(tome)
rownames(total) = gene_name

## read in more complete pheno data
pheno = read.csv(paste0(path, "sample_annotations.csv"),
    as.is = TRUE,
    row.names = 1)
pheno = DataFrame(pheno[colnames(total), ])

## convert
sce <-   SingleCellExperiment(list(counts = total), colData = pheno)

rowData(sce) = gene_map[rownames(sce), ]

## remove "exclude" cells
sce = sce[, sce$class_label != "Exclude"]
## save sce
save(sce, file = "rda/allen_snRNAseq_sce.Rdata")

## ## get pseudobulk
sce$BroadLayer = substr(sce$cortical_layer_label, 1, 2)
sce$SubjRegion = paste0(sce$region_label, "_", sce$external_donor_name_label)
sce$SubjRegionLayer = paste0(sce$SubjRegion, "_", sce$BroadLayer)

## broad class
sce$PseudoSample = sce$SubjRegionLayer


## split and combine
cIndexes = splitit(sce$PseudoSample)
umiComb <- sapply(cIndexes, function(ii)
    rowSums(assays(sce)$counts[, ii, drop = FALSE]))

phenoComb = colData(sce)[!duplicated(sce$PseudoSample),
    c(5:6, 11:21, 34:39)]

rownames(phenoComb) = phenoComb$PseudoSample
phenoComb = phenoComb[colnames(umiComb), ]
phenoComb = DataFrame(phenoComb)

sce_pseudobulk <-
    logNormCounts(SingleCellExperiment(
        list(counts = umiComb),
        colData = phenoComb,
        rowData = rowData(sce)
    ))

sce_pseudobulk$Layer_Class = paste0(sce_pseudobulk$cortical_layer_label,
    ":",
    sce_pseudobulk$class_label)

save(sce_pseudobulk, file = "rda/allen_scRNAseq_pseudobulked.Rdata")

###############################
##### get mean expression  ####

load("rda/allen_scRNAseq_pseudobulked.Rdata")
mat <- assays(sce_pseudobulk)$logcounts

## filter
gIndex = rowMeans(mat) > 1
mat_filter = mat[gIndex, ]

#####################
## Build a group model

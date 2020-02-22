###
#  module load conda_R/3.6.x
library(jaffelab)
library(Seurat)
library(scater)
library(DropletUtils)
library(limma)
library(lattice)
library(RColorBrewer)
library(Matrix)
library(parallel)

## read in data
pd = read.csv("mathys/snRNAseqPFC_BA10_biospecimen_metadata.csv", as.is=TRUE)
pheno = read.delim("mathys/filtered_column_metadata.txt", row.names = 1)
dat = readMM("mathys/filtered_count_matrix.mtx")
genes = read.delim("mathys/filtered_gene_row_names.txt",header=FALSE,as.is=TRUE)

## add names
rownames(dat) = genes$V1
colnames(dat) = rownames(pheno)

## get pseudobulk
pheno$individualID = pd$individualID[match(pheno$projid, pd$projid)]
pheno$PseudoSample = paste0(pheno$individualID, ":", pheno$Subcluster)

cIndexes = splitit(pheno$PseudoSample)
umiCombList <- mclapply(cIndexes, function(ii) {
	cat(".")
	rowSums(dat[, ii, drop = FALSE])
},mc.cores=12)
umiComb = do.call("cbind", umiCombList)

phenoComb = pheno[!duplicated(pheno$PseudoSample), 
	c("PseudoSample", "projid", "broad.cell.type",
		"Subcluster", "individualID")]
rownames(phenoComb) = phenoComb$PseudoSample
phenoComb = phenoComb[colnames(umiComb), ]
phenoComb = DataFrame(phenoComb)

## add more pheno
pd2 = read.csv("mathys/ROSMAP_Clinical_2019-05_v3.csv",as.is=TRUE)
pd2$age_death[pd2$age_death == "90+"] = 90
pd2$age_death = as.numeric(pd2$age_death)
pd2$Dx = factor(ifelse(pd2$age_first_ad_dx == "", "Control", "AD"),
	levels = c("Control", "AD"))
pd2 = pd2[match(phenoComb$individualID, pd2$individualID),]
phenoComb$Dx = pd2$Dx
phenoComb$age_death = pd2$age_death
phenoComb$msex = pd2$msex
phenoComb$race = pd2$race

sce_pseudobulk <-
    logNormCounts(SingleCellExperiment(
        list(counts = umiComb),
        colData = phenoComb,
        rowData = genes)
    )
save(sce_pseudobulk, file = "rda/mathys_pseudobulked.Rdata")

###############################
##### get mean expression  ####

load("rda/mathys_pseudobulked.Rdata")

mat_filter <- assays(sce_pseudobulk)$logcounts

#####################
## Build a group model
mod <- with(colData(sce_pseudobulk),
    model.matrix(~ 0 + Subcluster + Dx + age_death + msex + race))
colnames(mod) <- gsub('Subcluster', '', colnames(mod))

## get duplicate correlation
# corfit <- duplicateCorrelation(mat_filter, mod,
    # block = sce_pseudobulk$individualID)
# save(corfit, file = "rda/mathys_pseudobulked_dupCor.Rdata")
load("rda/mathys_pseudobulked_dupCor.Rdata")

## Next for each layer test that layer vs the rest
cell_idx <- splitit(sce_pseudobulk$Subcluster)

eb0_list_cell <- lapply(cell_idx, function(x) {
    res <- rep(0, ncol(sce_pseudobulk))
    res[x] <- 1
    m <- with(colData(sce_pseudobulk),
        model.matrix(~ res +
                Dx + age_death + msex + race))
    eBayes(
        lmFit(
            mat_filter,
            design = m,
            block = sce_pseudobulk$individualID,
            correlation = corfit$consensus.correlation
        )
    )
})
save(eb0_list_cell, file = "rda/mathys_pseudobulked_specific_Ts.Rdata")

##########
## Extract the p-values
load("rda/mathys_pseudobulked_specific_Ts.Rdata")

pvals0_contrasts_cell <- sapply(eb0_list_cell, function(x) {
    x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts_cell) = rownames(mat_filter)

t0_contrasts_cell <- sapply(eb0_list_cell, function(x) {
    x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts_cell) = rownames(mat_filter)
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
# Ast0   3061        1423        1131
# Ast1   2132        1015         838
# Ast2   1205         584         453
# Ast3    890         330         244
# End1    779         417         360
# End2    599         257         208
# Ex0    8870        2116        1323
# Ex1    9712        2477        1341
# Ex11   7266        1717        1052
# Ex12   2896         751         452
# Ex14   2469         669         444
# Ex2    1553         406         212
# Ex3    9417        2303        1238
# Ex4    4676        1016         545
# Ex5    8942        1832         981
# Ex6    3785         763         426
# Ex7    8828        1796        1041
# Ex8    3866        1883        1419
# Ex9    7081        1452         846
# In0    6866        1442         845
# In1    5615        1116         724
# In10    601         188         127
# In11    502         184         124
# In2    1345         437         304
# In3    1196         437         292
# In4     647         215         141
# In5    1129         402         283
# In6    1089         399         280
# In7    2057         514         342
# In8    1405         482         363
# In9    1035         393         291
# Mic0   2346        1271        1109
# Mic1   1727         975         831
# Mic2    693         416         374
# Mic3    563         245         188
# Oli0   6602        2531        1849
# Oli1   5674        2110        1615
# Oli3    679         285         204
# Oli4   1522         734         577
# Oli5   1306         644         521
# Opc0   3029        1206         920
# Opc1   2431        1048         815
# Opc2    435         135          90
# Per     947         474         387
############################
### correlate to layer?? ###
############################

###################
## load modeling outputs
load("rda/eb_contrasts.Rdata")
load("rda/eb0_list.Rdata")
load("rda/sce_layer.Rdata", verbose = TRUE)

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

mm = match(rowData(sce_layer)$gene_name, rownames(pvals0_contrasts_cell))

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
plot(hc)
cor_t_layer_toPlot = cor_t_layer[hc$order, c(1, 7:2)]

pdf("pdf/mathys_snRNAseq_overlap_heatmap.pdf", width = 11)
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

####################################
## relabel into spatial clusters ###
####################################

sce_pseudobulk$cellLayer = as.character(sce_pseudobulk$broad.cell.type)
sce_pseudobulk$cellLayer[sce_pseudobulk$Subcluster %in% paste0("Ex", c(0,2,4,6))] = "Ex_2/3"
sce_pseudobulk$cellLayer[sce_pseudobulk$Subcluster %in% paste0("Ex", c(1,5,14))] = "Ex_5"
sce_pseudobulk$cellLayer[sce_pseudobulk$Subcluster %in% paste0("Ex", c(7,8))] = "Ex_4"
sce_pseudobulk$cellLayer[sce_pseudobulk$Subcluster %in% paste0("Ex", c(3,11,12,9))] = "Ex_6"
sce_pseudobulk$cellLayer[sce_pseudobulk$Subcluster %in% paste0("In", c(0,7,9,11,2))] = "In_4/5"
sce_pseudobulk$cellLayer[sce_pseudobulk$Subcluster %in% paste0("In", c(10,3,6,1,4,5,8))] = "In_2/3"

## re-pseudo
sce_pseudobulk$PseudoSample_Layer = paste0(sce_pseudobulk$individualID, ":", sce_pseudobulk$cellLayer)

cIndexes_layer = splitit(sce_pseudobulk$PseudoSample_Layer)
umiComb <- sapply(cIndexes_layer, function(ii) {
	rowSums(assays(sce_pseudobulk)$counts[, ii, drop = FALSE])
})


phenoComb_layer = colData(sce_pseudobulk)
phenoComb_layer = phenoComb_layer[!duplicated(phenoComb_layer$PseudoSample_Layer),]
rownames(phenoComb_layer) = phenoComb_layer$PseudoSample_Layer
phenoComb_layer = phenoComb_layer[colnames(umiComb), ]
phenoComb_layer = DataFrame(phenoComb_layer)


sce_pseudobulk_layered <-
    logNormCounts(SingleCellExperiment(
        list(counts = umiComb),
        colData = phenoComb_layer,
        rowData = rowData(sce_pseudobulk)
    ))
save(sce_pseudobulk_layered, file = "rda/mathys_pseudobulked_layered.Rdata")

#####################
## Build a group model
mat_filter <- assays(sce_pseudobulk_layered)$logcounts
mod <- with(colData(sce_pseudobulk_layered),
    model.matrix(~ cellLayer * Dx + age_death + msex + race))
colnames(mod) <- gsub('cellLayer', '', colnames(mod))

# ## get duplicate correlation
corfit <- duplicateCorrelation(mat_filter, mod,
    block = sce_pseudobulk_layered$individualID)
save(corfit, file = "rda/mathys_pseudobulked_layered_dupCor.Rdata")
# load("rda/mathys_pseudobulked_layered_dupCor.Rdata")

fit = lmFit(mat_filter, design = mod,
            block = sce_pseudobulk_layered$individualID,
            correlation = corfit$consensus.correlation)
eb = eBayes(fit)

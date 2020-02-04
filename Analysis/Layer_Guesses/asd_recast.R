###
#  module load conda_R/3.6.x
library(jaffelab)
library(Seurat)
library(scater)
library(DropletUtils)
library(limma)

## read in data
pheno = read.delim("velmeshev/meta.tsv",row.names=1)
dat = read10xCounts("velmeshev/")
pheno = pheno[dat$Barcode,]
colData(dat) = DataFrame(pheno)

## get pseudobulk
dat$PseudoSample = paste0(dat$sample, ":", dat$cluster)

cIndexes = splitit(dat$PseudoSample)
umiComb <- sapply(cIndexes, function(ii)
        rowSums(assays(dat)$counts[, ii, drop = FALSE]))
		
phenoComb = colData(dat)[!duplicated(dat$PseudoSample),c(1:11,16)]
rownames(phenoComb) = phenoComb$PseudoSample
phenoComb = phenoComb[colnames(umiComb),]
phenoComb = DataFrame(phenoComb)

sce_pseudobulk <-
    logNormCounts(SingleCellExperiment(
        list(counts = umiComb),
        colData = phenoComb,
        rowData = rowData(dat)
    ))
save(sce_pseudobulk, file = "rda/velmeshev_pseudobulked.Rdata")

###############################
##### get mean expression  ####

load("rda/velmeshev_pseudobulked.Rdata")

mat <- assays(sce_pseudobulk)$logcounts

## filter 
gIndex = rowMeans(mat) > 0.1
mat_filter = mat[gIndex,]

#####################
## Build a group model
mod <- with(colData(sce_pseudobulk), 
		model.matrix( ~ 0 + cluster + region + age + sex + diagnosis))
colnames(mod) <- gsub('cluster', '', colnames(mod))

## get duplicate correlation
corfit <- duplicateCorrelation(mat_filter, mod, 
	block = sce_pseudobulk$individual)
save(corfit, file = "rda/velmeshev_pseudobulked_dupCor.Rdata")

## Next for each layer test that layer vs the rest
cell_idx <- splitit(sce_pseudobulk$cluster)

eb0_list_cell <- lapply(cell_idx, function(x) {
    res <- rep(0, ncol(sce_pseudobulk))
    res[x] <- 1
    m <- with(colData(sce_pseudobulk), 
		model.matrix( ~ res + 
			region + age + sex + diagnosis))
    eBayes(
        lmFit(
            mat_filter,
            design = m,
            block = sce_pseudobulk$individual,
            correlation = corfit$consensus.correlation
        )
    )
})
save(eb0_list_cell, file = "rda/velmeshev_pseudobulked_specific_Ts.Rdata")

## Extract the p-values
load("rda/velmeshev_pseudobulked_specific_Ts.Rdata")

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
    'FDRsig' = colSums(fdrs0_contrasts_cell< 0.05 & t0_contrasts_cell > 0),
    'Pval10-6sig' = colSums(pvals0_contrasts_cell < 1e-6 & t0_contrasts_cell > 0),
    'Pval10-8sig' = colSums(pvals0_contrasts_cell < 1e-8 & t0_contrasts_cell > 0)
)

                 # FDRsig Pval10.6sig Pval10.8sig
# AST-FB             3948        1800        1413
# AST-PP             5048        2272        1760
# Endothelial        3188        1781        1529
# IN-PV              7488        1422         854
# IN-SST             3659         851         583
# IN-SV2C            4668         988         660
# IN-VIP             7614        1463         946
# L2/3              10950        2435        1374
# L4                 7348        1266         755
# L5/6               5369        1301         816
# L5/6-CC            7877        1300         688
# Microglia          2200        1250        1068
# Neu-mat            1896         397         218
# Neu-NRGN-I         3376        1302         957
# Neu-NRGN-II        4138        1922        1421
# Oligodendrocytes   4187        1893        1515
# OPC                5812        2001        1468
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

mm = match(rownames(pvals0_contrasts), rownames(pvals0_contrasts_cell))

pvals0_contrasts = pvals0_contrasts[!is.na(mm),]
t0_contrasts = t0_contrasts[!is.na(mm),]
fdrs0_contrasts = fdrs0_contrasts[!is.na(mm),]

pvals0_contrasts_cell = pvals0_contrasts_cell[mm[!is.na(mm)],]
t0_contrasts_cell = t0_contrasts_cell[mm[!is.na(mm)],]
fdrs0_contrasts_cell = fdrs0_contrasts_cell[mm[!is.na(mm)],]

cor_t = cor(t0_contrasts_cell,t0_contrasts)
signif(cor_t,2)

### just layer specific genes from ones left
layer_specific_indices = mapply(function(t,p) {
	oo = order(t, decreasing=TRUE)[1:100]
	}, as.data.frame(t0_contrasts), as.data.frame(pvals0_contrasts))
layer_ind = unique(as.numeric(layer_specific_indices))

cor_t_layer = cor(t0_contrasts_cell[layer_ind,],
	t0_contrasts[layer_ind,])
signif(cor_t_layer,3)

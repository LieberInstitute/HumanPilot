##
library(limma)
library(jaffelab)

## load outputs
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

## finding layer specific genes
layer_specific_indices = mapply(function(t,p) {
	oo = order(t, decreasing=TRUE)[1:25]
	}, as.data.frame(t0_contrasts), as.data.frame(pvals0_contrasts))
layer_ind = as.numeric(layer_specific_indices)

t0_contrasts_sig = t0_contrasts[layer_ind,]
fdrs0_contrasts_sig = fdrs0_contrasts[layer_ind,]
pvals0_contrasts_sig = pvals0_contrasts[layer_ind,]

## load layer level
load("rda/sce_layer.Rdata")

## means by layer
layer_indices = splitit(sce_layer$layer_guess)
layer_log = assays(sce_layer)$logcounts
layer_means = sapply(layer_indices, function(ii) rowMeans(layer_log[,ii]))
layer_means_sig = layer_means[rownames(pvals0_contrasts_sig),]
layer_log_sig = layer_log[rownames(pvals0_contrasts_sig),]
image(layer_log_sig)
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
load(here('Analysis', 'Layer_Guesses',
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

## filter to same genes
sce = sce[rownames(t0_contrasts),]

#########################
### do correlations #####
#########################

## get out counts
sce_logcounts = assays(sce)$logcounts

## using all genes
cc_spots = cor(as.matrix(sce_logcounts),t0_contrasts)
cc_spots = as.data.frame(cc_spots)
cc_spots$subject = sce$subject
cc_spots$layer_guess = sce$layer_guess
cc_spots$cell_count = sce$cell_count

pdf("layer_correlations.pdf",w=12)
par(mar=c(10,6,2,2))
boxplot(WM ~subject*layer_guess,data=cc_spots, las=3,xlab="",ylab="Correlation")
boxplot(Layer1 ~subject*layer_guess,data=cc_spots, las=3,xlab="",ylab="Correlation")
boxplot(Layer2 ~subject*layer_guess,data=cc_spots, las=3,xlab="",ylab="Correlation")
boxplot(Layer3 ~subject*layer_guess,data=cc_spots, las=3,xlab="",ylab="Correlation")
boxplot(Layer4 ~subject*layer_guess,data=cc_spots, las=3,xlab="",ylab="Correlation")
boxplot(Layer5 ~subject*layer_guess,data=cc_spots, las=3,xlab="",ylab="Correlation")
boxplot(Layer6 ~subject*layer_guess,data=cc_spots, las=3,xlab="",ylab="Correlation")

###
## using divergent genes
g = rownames(t0_contrasts_sig)
cc_spots_sig = cor(as.matrix(sce_logcounts[g,]),t0_contrasts[g,])
cc_spots_sig = as.data.frame(cc_spots_sig)
cc_spots_sig$subject = sce$subject
cc_spots_sig$layer_guess = sce$layer_guess
cc_spots_sig$cell_count = sce$cell_count

boxplot(WM ~subject*layer_guess,data=cc_spots_sig, las=3,xlab="",ylab="Correlation")
boxplot(Layer1 ~subject*layer_guess,data=cc_spots_sig, las=3,xlab="",ylab="Correlation")
boxplot(Layer2 ~subject*layer_guess,data=cc_spots_sig, las=3,xlab="",ylab="Correlation")
boxplot(Layer3 ~subject*layer_guess,data=cc_spots_sig, las=3,xlab="",ylab="Correlation")
boxplot(Layer4 ~subject*layer_guess,data=cc_spots_sig, las=3,xlab="",ylab="Correlation")
boxplot(Layer5 ~subject*layer_guess,data=cc_spots_sig, las=3,xlab="",ylab="Correlation")
boxplot(Layer6 ~subject*layer_guess,data=cc_spots_sig, las=3,xlab="",ylab="Correlation")

## one cell
boxplot(WM ~subject*layer_guess,data=cc_spots_sig, subset = cell_count == 1, las=3,xlab="",ylab="Correlation")
boxplot(Layer1 ~subject*layer_guess,data=cc_spots_sig, subset = cell_count == 1,las=3,xlab="",ylab="Correlation")
boxplot(Layer2 ~subject*layer_guess,data=cc_spots_sig, subset = cell_count == 1,las=3,xlab="",ylab="Correlation")
boxplot(Layer3 ~subject*layer_guess,data=cc_spots_sig, subset = cell_count == 1,las=3,xlab="",ylab="Correlation")
boxplot(Layer4 ~subject*layer_guess,data=cc_spots_sig, subset = cell_count == 1,las=3,xlab="",ylab="Correlation")
boxplot(Layer5 ~subject*layer_guess,data=cc_spots_sig, subset = cell_count == 1,las=3,xlab="",ylab="Correlation")
boxplot(Layer6 ~subject*layer_guess,data=cc_spots_sig, subset = cell_count == 1,las=3,xlab="",ylab="Correlation")

boxplot(WM ~ cell_count,data=cc_spots_sig)
boxplot(Layer1 ~ cell_count,data=cc_spots_sig)
boxplot(Layer2 ~ cell_count,data=cc_spots_sig)
boxplot(Layer3 ~ cell_count,data=cc_spots_sig)
boxplot(Layer4 ~ cell_count,data=cc_spots_sig)
boxplot(Layer5 ~ cell_count,data=cc_spots_sig)
boxplot(Layer6 ~ cell_count,data=cc_spots_sig)
dev.off()

save(cc_spots, cc_spots_sig, file = "rda/correlations_of_layer_stats_to_spots.Rdata")

### try spot level


sce_image_grid_gene(
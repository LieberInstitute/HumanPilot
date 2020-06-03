##
library('SingleCellExperiment')
library('here')
library('readxl')
library('Polychrome')
library('rafalib')
library('sessioninfo')
library('WGCNA')
library('lmerTest')

## multithread
allowWGCNAThreads(6)

## Functions derived from this script, to make it easier to resume the work
sce_layer_file <-
    here('Analysis', 'Layer_Guesses', 'rda', 'sce_layer.Rdata')
if (file.exists(sce_layer_file))
    load(sce_layer_file, verbose = TRUE)
source(here('Analysis', 'Layer_Guesses', 'layer_specificity_functions.R'),
    echo = TRUE)
	
	
#########################
## get power
geneExprs = assays(sce_layer)$logcounts

powers <- c(1:10, seq(from = 12, to=20, by=2))
sftthresh <- pickSoftThreshold(t(geneExprs), powerVector = powers,
                               networkType = "signed", verbose = 5)
cat(sftthresh$powerEstimate) # 7

## run wgcna
net = blockwiseModules(t(geneExprs), power = sftthresh$powerEstimate,
					networkType = "signed", minModuleSize = 30,
					corType="bicor", reassignThreshold = 0, 
					mergeCutHeight = 0.25, numericLabels = TRUE, 
					pamRespectsDendro = FALSE,verbose = 5)
save(net, sftthresh, file = "rda/constructed_network_signed_bicor_layerLevel.rda")


### checks
load( "rda/constructed_network_signed_bicor_layerLevel.rda")

table(net$colors)
    # 0     1     2     3     4     5     6     7     8     9    10    11    12
# 10789  2476  1621  1360  1047   932   880   744   414   373   349   344   325
   # 13    14    15    16    17    18
  # 190   174   163    60    46    44

sce_layer$layer = gsub("ayer", "", as.character(sce_layer$layer_guess))
sce_layer$layer= factor(sce_layer$layer, 
	levels = c(paste0("L", 1:6), "WM"))

MEs = net$MEs
MEs = MEs[,paste0("ME",names(table(net$colors)))]

# boxplots
mypar(3,7)
for(i in 1:ncol(MEs)) {
	boxplot(MEs[,i] ~ sce_layer$layer,ylab = "",
		xlab = "", las=3, main = colnames(MEs)[i])
}

## mean expression per module
means = rowMeans(geneExprs)
signif(tapply(means, net$colors, mean),2)

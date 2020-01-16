library('SingleCellExperiment')
library('here')
library('jaffelab')
library('scater')
library('scran')
library('pheatmap')
library('sessioninfo')

## Load data
load(here(
    'Analysis',
    'Human_DLPFC_Visium_processedData_sce_scran.Rdata'
))

## Load layer guesses
load(here(
    'Analysis', 'Layer_Guesses',
    'layer_guess_tab.Rdata'
))

## Add layer guesses to the sce object
sce$layer_guess <- NA
m <- match(sce$key, layer_guess_tab$key)
table(is.na(m))
# FALSE  TRUE
# 47329   352
sce$layer_guess[!is.na(m)] <- layer_guess_tab$layer[m[!is.na(m)]]

## Check layer guesses across images
options(width = 120)
with(colData(sce), addmargins(table(layer_guess, sample_name, useNA = 'ifany')))
#            sample_name
# layer_guess 151507 151508 151509 151510 151669 151670 151671 151672 151673 151674 151675 151676   Sum
#   Layer 1      817    866   1189   1180      0      0      0      0    273    380    328    289  5322
#   Layer 2      305    295    602    650      0      0      0      0    253    224    275    254  2858
#   Layer 2/3      0      0      0      0   2141   2175   1918   1575      0      0      0      0  7809
#   Layer 3     1215   1385   1884   1774      0      0      0      0    989    924    771    836  9778
#   Layer 4      369    373    369    318    364    211    245    304    218    247    275    254  3547
#   Layer 5      675    737    363    310    510    581    721    728    673    621    732    649  7300
#   Layer 6      486    525    215    179    391    308    760    882    692    614    533    616  6201
#   WM           354    200    166    184    230    209    449    399    513    625    652    533  4514
#   <NA>           5      3      1     39     25     14     17    127     28     38     26     29   352
#   Sum         4226   4384   4789   4634   3661   3498   4110   4015   3639   3673   3592   3460 47681

## Drop the layer guess NAs for now
sce <- sce[, !is.na(sce$layer_guess)]
dim(sce)
# [1] 33538 47329

## Next, re-label "Layer 2/3" as "Layer 3" for now
## (there's more layer 3 in the other samples than 2 anyway)
sce$layer_guess[sce$layer_guess == 'Layer 2/3'] <- 'Layer 3'

## Make it into a factor with WM as the reference
sce$layer_guess <- factor(sce$layer_guess, levels = c('WM', paste('Layer', 1:6)))

## Check again
with(colData(sce), addmargins(table(layer_guess, sample_name, useNA = 'ifany')))
#            sample_name
# layer_guess 151507 151508 151509 151510 151669 151670 151671 151672 151673 151674 151675 151676   Sum
#     WM         354    200    166    184    230    209    449    399    513    625    652    533  4514
#     Layer 1    817    866   1189   1180      0      0      0      0    273    380    328    289  5322
#     Layer 2    305    295    602    650      0      0      0      0    253    224    275    254  2858
#     Layer 3   1215   1385   1884   1774   2141   2175   1918   1575    989    924    771    836 17587
#     Layer 4    369    373    369    318    364    211    245    304    218    247    275    254  3547
#     Layer 5    675    737    363    310    510    581    721    728    673    621    732    649  7300
#     Layer 6    486    525    215    179    391    308    760    882    692    614    533    616  6201
#     Sum       4221   4381   4788   4595   3636   3484   4093   3888   3611   3635   3566   3431 47329
    
## Split spots by layer and image
## Based on here('Analysis', 'sce_scran.R')
## /dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/sce_scran.R
layerIndexes <- splitit(paste0(sce$sample_name, '_', sce$layer_guess))

## Collapse UMIs
umiComb <-
    sapply(layerIndexes, function(ii)
        rowSums(assays(sce)$counts[, ii, drop = FALSE]))
dim(umiComb)
# [1] 33538    76
## That's because 2 layers are absent in subject 2
stopifnot(12 * 7 - 2 * 4 == ncol(umiComb))

## Get a sample-specific size factors, instead of sample/layer size factors
umiComb_sample <-
    sapply(splitit(ss(colnames(umiComb), '_', 1)), function(ii) {
        rowSums(umiComb[, ii, drop = FALSE])
    })
umiComb_sample_size_fac <- librarySizeFactors(umiComb_sample)

umiComb_sample_size_fac_layer <-
    rep(umiComb_sample_size_fac, lengths(splitit(ss(
        colnames(umiComb), '_', 1
    ))))
names(umiComb_sample_size_fac_layer) <- colnames(umiComb)

## Build a data.frame with the pheno data by layer
layer_df <- data.frame(
    sample_name = ss(colnames(umiComb), '_', 1),
    layer_guess = factor(ss(colnames(umiComb), '_', 2), levels = levels(sce$layer_guess)),
    stringsAsFactors = FALSE
)
m_layer <- match(layer_df$sample_name, sce$sample_name)
layer_df$subject <- sce$subject[m_layer]
layer_df$position <- sce$position[m_layer]
layer_df$replicate <- sce$replicate[m_layer]
layer_df$subject_position <- sce$subject_position[m_layer]
rownames(layer_df) <- colnames(umiComb)

## Build a new sce object with the library size factors
## by sample instead of sample/layer
sce_layer <-
    logNormCounts(
        SingleCellExperiment(list(counts = umiComb),
        colData = layer_df,
        rowData = rowData(sce)),
        size_factors = umiComb_sample_size_fac_layer
    )


## Try using the scran::findMarkers() function instead of coding the
## looping code myself
## More at https://osca.bioconductor.org/marker-detection.html

## 
markers_layer_any <- findMarkers(sce_layer, sce_layer$layer_guess, pval.type = 'all', direction = 'up', block = sce_layer$subject_position)


## I downloaded the current git version (026edd8eb68fa0479769450473a31795ae50c742)
## to obtain the code for this function
## git clone https://git.bioconductor.org/packages/scran
getMarkerEffects <- function(x, prefix="logFC", strip=TRUE) {
    regex <- paste0("^", prefix, "\\.")
    i <- grep(regex, colnames(x))
    out <- as.matrix(x[,i])

    if (strip) {
        colnames(out) <- sub(regex, "", colnames(out))
    }
    out
}


to_symbols <- function(x) {
    m <- match(rownames(x), rownames(rowData(sce_layer)))
    rownames(x) <- rowData(sce_layer)$gene_name[m]
    return(x)
}



pdf('test2.pdf', height = 20, useDingbats = FALSE)
lapply(seq_along(markers_layer_any), function(chosen) {
    interesting <- to_symbols(markers_layer_any[[chosen]])
    # best.set <- interesting[interesting$Top <= 6,] ## for pval.type = 'any'
    best.set <- head(interesting, 30) ## for pval.type == 'all'
    logFCs <- getMarkerEffects(best.set)
    print(pheatmap(logFCs, breaks=seq(-5, 5, length.out=101), main = names(markers_layer_any)[chosen]))
    return(NULL)
})
dev.off()


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
library('SingleCellExperiment')
library('here')
library('readxl')
library('limma')
library('sessioninfo')

dir.create('pdf', showWarnings = FALSE)
dir.create('rda', showWarnings = FALSE)

## Load data
load(here(
    'Analysis',
    'Human_DLPFC_Visium_processedData_sce_scran.Rdata'
))

## Functions derived from this script, to make it easier to resume the work
sce_layer_file <-
    here('Analysis', 'Layer_Guesses', 'rda', 'sce_layer.Rdata')
if (file.exists(sce_layer_file))
    load(sce_layer_file, verbose = TRUE)
source(here('Analysis', 'Layer_Guesses', 'layer_specificity_functions.R'))

## For plotting
source(here('Analysis', 'spatialLIBD_global_plot_code.R'))
genes <- paste0(rowData(sce)$gene_name, '; ', rowData(sce)$gene_id)


## mean Xk unique molecular indices (UMIs) and mean Xk genes per spot
summary(sce$sum_umi)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   17    2035    3034    3462    4407   20600

summary(sce$sum_gene)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#   16    1178    1631    1734    2176    6035


## chrM genes
ix_mito <- grep("^MT-", rowData(sce)$gene_name)
rowData(sce)$gene_name[ix_mito]
# [1] "MT-ND1"  "MT-ND2"  "MT-CO1"  "MT-CO2"  "MT-ATP8" "MT-ATP6" "MT-CO3"
# [8] "MT-ND3"  "MT-ND4L" "MT-ND4"  "MT-ND5"  "MT-ND6"  "MT-CYB"

## Should save this on the sce object later
expr_total <- colSums(assays(sce)$counts)
## Actually, we already had this
identical(sce$sum_umi, expr_total)
# [1] TRUE
expr_chrM <- colSums(assays(sce)$counts[ix_mito,])
expr_chrM_ratio <- expr_chrM / expr_total
summary(expr_chrM_ratio)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.04853 0.15465 0.18442 0.18554 0.21521 0.44156


## Visualize this data at the spot level
## In the future we could customize the colors if we want to
sce$expr_total <- expr_total
sce$expr_chrM <- expr_chrM
sce$expr_chrM_ratio <- expr_chrM_ratio
sce_image_grid_gene(
    sce,
    geneid = 'expr_total',
    spatial = TRUE,
    minCount = 0,
    pdf_file = 'pdf/spot_expr_total.pdf'
)
sce_image_grid_gene(
    sce,
    geneid = 'expr_chrM',
    spatial = TRUE,
    minCount = 0,
    pdf_file = 'pdf/spot_expr_chrM.pdf'
)
sce_image_grid_gene(
    sce,
    geneid = 'expr_chrM_ratio',
    spatial = TRUE,
    minCount = 0,
    pdf_file = 'pdf/spot_expr_chrM_ratio.pdf'
)


## Repeat at the layer-level
# ix_mito_layer <- grep("^MT-", rowData(sce_layer)$gene_name)
# expr_total_layer <- colSums(assays(sce_layer)$counts)
# expr_chrM_layer <-
#     colSums(assays(sce_layer)$counts[ix_mito_layer,])
# expr_chrM_ratio_layer <- expr_chrM_layer / expr_total_layer
# summary(expr_chrM_ratio_layer)
## Err, it's all 0 because we already dropped chrM by this point :P

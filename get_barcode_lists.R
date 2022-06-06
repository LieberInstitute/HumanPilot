####
###
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
library('reshape2')
library('lmerTest')

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
load(here('Analysis', 'Layer_Guesses', 'rda',
    'layer_guess_tab.Rdata'))

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
sce_original <- sce
sce <- sce[, !is.na(sce$layer_guess)]
dim(sce)
sce$layer_guess[sce$layer_guess == 'Layer 2/3'] <- 'Layer 3'


## make big table
pd = colData(sce)[,c("barcode", "sample_name", "layer_guess")]
pd$layer_guess = gsub("ayer ", "", pd$layer_guess)
pd$barcode = as.character(pd$barcode)

sample_tab = pd[!duplicated(pd[,-1]),-1]
sample_tab$BAM = paste0("/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/10X/",
	sample_tab$sample_name, "/", sample_tab$sample_name, "_mRNA.bam")

write.table(as.data.frame(pd), "10X/barcode_level_layer_map.tsv",
	row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
write.table(as.data.frame(sample_tab), "10X/sample_level_layer_map.tsv",
	row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
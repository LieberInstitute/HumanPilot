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

dir.create('pdf', showWarnings = FALSE)
dir.create('rda', showWarnings = FALSE)

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
# [1] 33538 47329

## Next, re-label "Layer 2/3" as "Layer 3" for now
## (there's more layer 3 in the other samples than 2 anyway)
sce$layer_guess[sce$layer_guess == 'Layer 2/3'] <- 'Layer 3'

## Make it into a factor with WM as the reference
## and remove spaces
sce$layer_guess <-
    factor(gsub(' ', '', sce$layer_guess), levels = c('WM', paste0('Layer', 1:6)))

## Check again
with(colData(sce), addmargins(table(layer_guess, sample_name, useNA = 'ifany')))
#            sample_name
# layer_guess 151507 151508 151509 151510 151669 151670 151671 151672 151673 151674 151675 151676   Sum
#     WM         354    200    166    184    230    209    449    399    513    625    652    533  4514
#     Layer1    817    866   1189   1180      0      0      0      0    273    380    328    289  5322
#     Layer2    305    295    602    650      0      0      0      0    253    224    275    254  2858
#     Layer3   1215   1385   1884   1774   2141   2175   1918   1575    989    924    771    836 17587
#     Layer4    369    373    369    318    364    211    245    304    218    247    275    254  3547
#     Layer5    675    737    363    310    510    581    721    728    673    621    732    649  7300
#     Layer6    486    525    215    179    391    308    760    882    692    614    533    616  6201
#     Sum       4221   4381   4788   4595   3636   3484   4093   3888   3611   3635   3566   3431 47329

## Split spots by layer and image
## Based on here('Analysis', 'sce_scran.R')
## /dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/sce_scran.R
layerIndexes <-
    splitit(paste0(sce$sample_name, '_', sce$layer_guess))

## Collapse UMIs
umiComb <-
    sapply(layerIndexes, function(ii)
        rowSums(assays(sce)$counts[, ii, drop = FALSE]))
dim(umiComb)
# [1] 33538    76
## That's because 2 layers are absent in subject 2
stopifnot(12 * 7 - 2 * 4 == ncol(umiComb))

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

## Compare the size factors
pdf('pdf/size_factors_comparison.pdf', useDingbats = FALSE)
plot(
    librarySizeFactors(umiComb),
    umiComb_sample_size_fac_layer,
    bg = Polychrome::palette36.colors(7)[as.integer(layer_df$layer_guess)],
    pch = 21,
    xlab = 'Default size factors',
    ylab = 'Sample size factors (repeated across layers)',
    main = 'Compare size factors'
)
legend(
    x = 1,
    y = 1.5,
    legend = levels(layer_df$layer_guess),
    col =  Polychrome::palette36.colors(7),
    lwd = 3,
    bty = 'n',
    ncol = 2
)
dev.off()

## There are layer differences, so it's likely best to use the default
## size factors instead

## Build a new sce object
sce_layer <-
    logNormCounts(SingleCellExperiment(
        list(counts = umiComb),
        colData = layer_df,
        rowData = rowData(sce)
    ))
# ,size_factors = umiComb_sample_size_fac_layer)

## From Lukas's code
## https://github.com/LieberInstitute/HumanPilot/commit/0bd87cb5cd69863ceafc343f5456dcafca0dccb6#diff-6c1a4dea06c53f935341cef1d3bccbfdR158
## for identifying the mitochondrial genes
ix_mito <- grep("^MT-", rowData(sce_layer)$gene_name)

## Several of the mitochondrial genes are among the most expressed ones
pdf('pdf/highest_expressed.pdf', useDingbats = FALSE)
plotHighestExprs(
    sce_layer,
    exprs_values = "counts",
    colour_cells_by = 'layer_guess',
    feature_names_to_plot = 'gene_name',
) + ggplot2::ggtitle('All genes') + ggplot2::xlab('Percent among counts') +
    ggplot2::scale_colour_manual(values =  unname(Polychrome::palette36.colors(7)),
        name = 'Layer')
plotHighestExprs(
    sce_layer,
    exprs_values = "counts",
    colour_cells_by = 'layer_guess',
    feature_names_to_plot = 'gene_name',
    drop_features = seq_len(nrow(sce_layer))[-ix_mito]
) + ggplot2::ggtitle('Only mitochondrial genes') +
    ggplot2::xlab('Percent among counts') +
    ggplot2::scale_colour_manual(values =  unname(Polychrome::palette36.colors(7)),
        name = 'Layer')
plotHighestExprs(
    sce_layer,
    exprs_values = "counts",
    colour_cells_by = 'layer_guess',
    feature_names_to_plot = 'gene_name',
    drop_features = ix_mito,
) + ggplot2::ggtitle('Without mitochondrial genes') +
    ggplot2::xlab('Percent among counts') +
    ggplot2::scale_colour_manual(values =  unname(Polychrome::palette36.colors(7)),
        name = 'Layer')
dev.off()


pdf('pdf/mt_expression.pdf', useDingbats = FALSE)
## Hack my way around it so I can show the gene symbols
x <-
    plotExpression(sce_layer, rownames(sce_layer)[ix_mito], x = 'layer_guess')
x$data$Feature <-
    rowData(sce_layer)$gene_name[match(as.character(x$data$Feature), rownames(sce_layer))]
print(x)

## Different yet similar hack
x <-
    plotExpression(sce_layer, rownames(sce_layer)[ix_mito], colour_by = 'layer_guess')
x$data$X <-
    rowData(sce_layer)$gene_name[match(as.character(x$data$X), rownames(sce_layer))]
print(x +
        ggplot2::scale_fill_manual(
            values =  unname(Polychrome::palette36.colors(7)),
            name = 'Layer'
        ))
dev.off()
rm(x)

## Drop mitochondrial genes
sce_layer <- sce_layer[-ix_mito, ]


## Find which genes to drop due to low expression values
sce_layer_avg <- calculateAverage(sce_layer)
summary(sce_layer_avg)
#  Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.000    0.007    1.186   53.215   31.427 9507.091

## Results without dropping the mitochondrial genes:
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.00     0.01     1.19    64.31    31.50 59084.15


sce_layer_avg_logcounts <-
    calculateAverage(sce_layer, exprs_values = 'logcounts')
summary(sce_layer_avg_logcounts)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.000000  0.003973  1.379983  4.770529  9.678251 26.245639

## Results without dropping the mitochondrial genes:
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.000000  0.003973  1.384647  4.779067  9.683280 30.899169

# hist(sce_layer_avg)
# hist(sce_layer_avg_logcounts, breaks = 50)

addmargins(
    table(
        'avg log > 0.1' = sce_layer_avg_logcounts > 0.1,
        'avg log > 0.05' = sce_layer_avg_logcounts > 0.05,
        'avg log > 0' = sce_layer_avg_logcounts > 0
    )
)
# , , avg log > 0 = FALSE
#
#              avg log > 0.05
# avg log > 0.1 FALSE  TRUE   Sum
#         FALSE  7943     0  7943
#         TRUE      0     0     0
#         Sum    7943     0  7943
#
# , , avg log > 0 = TRUE
#
#              avg log > 0.05
# avg log > 0.1 FALSE  TRUE   Sum
#         FALSE  2746  1015  3761
#         TRUE      0 21821 21821
#         Sum    2746 22836 25582
#
# , , avg log > 0 = Sum
#
#              avg log > 0.05
# avg log > 0.1 FALSE  TRUE   Sum
#         FALSE 10689  1015 11704
#         TRUE      0 21821 21821
#         Sum   10689 22836 33525

## From scater's vignette at
## http://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/overview.html#333_subsetting_by_row
keep_feature <- nexprs(sce_layer, byrow = TRUE) > 0

## The avg log counts cutoff would be more strict
addmargins(table(keep_feature, 'avg log > 0.05' = sce_layer_avg_logcounts > 0.05))
#             avg log > 0.05
# keep_feature FALSE  TRUE   Sum
#        FALSE  7943     0  7943
#        TRUE   2746 22836 25582
#        Sum   10689 22836 33525

## From scater's vignette:
# detected: the percentage of cells with non-zero counts for each gene.
per.feat <- perFeatureQCMetrics(sce_layer)
summary(per.feat$detected)
#  Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.000   1.316  55.263  51.674 100.000 100.000

addmargins(
    table(
        'Detected >5%' = per.feat$detected > 5,
        'avg log > 0.05' = sce_layer_avg_logcounts > 0.05,
        keep_feature
    )
)
# , , keep_feature = FALSE
#
#             avg log > 0.05
# Detected >5% FALSE  TRUE   Sum
#        FALSE  7943     0  7943
#        TRUE      0     0     0
#        Sum    7943     0  7943
#
# , , keep_feature = TRUE
#
#             avg log > 0.05
# Detected >5% FALSE  TRUE   Sum
#        FALSE  2389   505  2894
#        TRUE    357 22331 22688
#        Sum    2746 22836 25582
#
# , , keep_feature = Sum
#
#             avg log > 0.05
# Detected >5% FALSE  TRUE   Sum
#        FALSE 10332   505 10837
#        TRUE    357 22331 22688
#        Sum   10689 22836 33525

pdf('pdf/avg_logcounts_VS_detected_percent.pdf', useDingbats = FALSE)
plot(
    sce_layer_avg_logcounts,
    log10(per.feat$detected),
    bg = add.alpha('black', 1 / 10),
    pch = 21,
    col = NA,
    xlab = 'Average logcounts',
    ylab = 'Layer/image % detected (log10)'
)
abline(v = 0.05,
    col = 'red',
    lwd = 2,
    lty = 2)
abline(
    h = log10(5),
    col = 'blue',
    lwd = 2,
    lty = 2
)

plot(
    sce_layer_avg_logcounts,
    per.feat$detected,
    xlim = c(0, 0.1),
    ylim = c(0, 10),
    bg = add.alpha('black', 1 / 10),
    pch = 21,
    col = NA,
    xlab = 'Average logcounts',
    ylab = 'Layer/image % detected'
)
abline(v = 0.05,
    col = 'red',
    lwd = 2,
    lty = 2)
abline(h = 5,
    col = 'blue',
    lwd = 2,
    lty = 2)

plot(
    sce_layer_avg_logcounts,
    log10(per.feat$detected),
    xlim = c(0.05, max(sce_layer_avg_logcounts)),
    ylim = log10(c(5, 100)),
    bg = add.alpha('black', 1 / 10),
    pch = 21,
    col = NA,
    xlab = 'Average logcounts',
    ylab = 'Layer/image % detected  (log10)'
)
abline(v = 0.05,
    col = 'red',
    lwd = 2,
    lty = 2)
abline(
    h = log10(5),
    col = 'blue',
    lwd = 2,
    lty = 2
)
dev.off()

## The 5% detected means that the gene must be expressed in at least 4/76 layer/image combinations
4 / 76
# [1] 0.05263158

## Drop genes
selected_genes <-
    which(per.feat$detected > 5 & sce_layer_avg_logcounts > 0.05)
length(selected_genes)
# [1] 22331
sce_layer <- sce_layer[selected_genes,]
dim(sce_layer)
# [1] 22331    76

## For comparison, try with the other library size factors
sce_layer_sfac <-
    logNormCounts(SingleCellExperiment(
        list(counts = umiComb),
        colData = layer_df,
        rowData = rowData(sce)
    ),
        size_factors = umiComb_sample_size_fac_layer)[-ix_mito, ][selected_genes,]

## From scater's vignette at
## http://bioconductor.org/packages/release/bioc/vignettes/scater/inst/doc/overview.html#34_variable-level_qc
vars <- getVarianceExplained(
    sce_layer,
    variables = c(
        'sample_name',
        'layer_guess',
        'subject',
        'subject_position',
        'position'
    )
)
vars_sfac <- getVarianceExplained(
    sce_layer_sfac,
    variables = c(
        'sample_name',
        'layer_guess',
        'subject',
        'subject_position',
        'position'
    )
)
pdf('pdf/gene_explanatory_vars.pdf', useDingbats = FALSE)
plotExplanatoryVariables(vars) + ggplot2::ggtitle('With default sample size factors; no MT')
## Hack it so the colors are the same in both plots
plotExplanatoryVariables(vars_sfac) +
    ggplot2::ggtitle('With sample size factors (repeated across layers); no MT') +
    ggplot2::scale_colour_manual(values =  unname(scater:::.get_palette('tableau10medium')[c(2, 1, 3:5)]),
        name = '')
dev.off()
## From the above plots it indeed looks like using the default size factors is best

## From here('Analysis', 'sce_scran.R')
## Find the highly variable genes
dec <- modelGeneVar(sce_layer,
    block = sce_layer$subject_position)

pdf('pdf/modelGeneVar.pdf', useDingbats = FALSE)
mapply(function(block, blockname) {
    plot(
        block$mean,
        block$total,
        xlab = "Mean log-expression",
        ylab = "Variance",
        main = blockname
    )
    curve(metadata(block)$trend(x),
        col = "blue",
        add = TRUE)
}, dec$per.block, names(dec$per.block))
dev.off()

top.hvgs <- getTopHVGs(dec, prop = 0.1)
length(top.hvgs)
# [1] 1280

## Run PCA with a reduced number of components to avoid errors
## (since we are starting with 76 columns anyway instead of thousands)
sce_layer <-
    runPCA(sce_layer, subset_row = top.hvgs, ncomponents = 20)

## Default perplexity is 15
mat <-
    scater:::.get_mat_from_sce(
        sce_layer,
        exprs_values = 'logcounts',
        dimred = 'PCA',
        n_dimred = NULL
    )
dim(mat)
min(50, floor(nrow(mat) / 5))
# [1] 15

## Run TSNE and UMAP
set.seed(20200121)
sce_layer <-
    runTSNE(sce_layer,
        dimred = 'PCA',
        name = 'TSNE_perplexity5',
        perplexity = 5)
set.seed(20200121)
sce_layer <-
    runTSNE(sce_layer,
        dimred = 'PCA',
        name = 'TSNE_perplexity15',
        perplexity = 15)
set.seed(20200121)
sce_layer <-
    runTSNE(sce_layer,
        dimred = 'PCA',
        name = 'TSNE_perplexity20',
        perplexity = 20)
set.seed(20200121)
sce_layer <-
    runUMAP(sce_layer, dimred = 'PCA', name = 'UMAP_neighbors15')

sce_layer$layer_guess_reordered <-
    factor(sce_layer$layer_guess, levels = c(paste0('Layer', 1:6), 'WM'))
pdf('pdf/reduced_dim_PCA.pdf', useDingbats = FALSE)
plotReducedDim(
    sce_layer,
    dimred = 'PCA',
    colour_by = 'layer_guess_reordered',
    theme_size = 20,
    point_size = 5
) +
    ggplot2::scale_fill_manual(values =  unname(Polychrome::palette36.colors(7))[c(2:7, 1)],
        name = 'Layer')
plotReducedDim(
    sce_layer,
    dimred = 'PCA',
    colour_by = 'sample_name',
    theme_size = 20,
    point_size = 5
)
plotReducedDim(
    sce_layer,
    dimred = 'PCA',
    colour_by = 'subject_position',
    theme_size = 20,
    point_size = 5
)
plotReducedDim(
    sce_layer,
    dimred = 'PCA',
    colour_by = 'subject',
    theme_size = 20,
    point_size = 5
)
plotReducedDim(
    sce_layer,
    dimred = 'PCA',
    colour_by = 'position',
    theme_size = 20,
    point_size = 5
)
dev.off()

## Just to compare the PCA on the data using library size factors
## by sample repeated across layers
sce_layer_sfac <-
    runPCA(sce_layer_sfac,
        subset_row = top.hvgs,
        ncomponents = 20)
pdf('pdf/reduced_dim_PCA_sfact.pdf', useDingbats = FALSE)
plotReducedDim(
    sce_layer_sfac,
    dimred = 'PCA',
    colour_by = 'layer_guess',
    theme_size = 20
) +
    ggplot2::scale_fill_manual(values =  unname(Polychrome::palette36.colors(7)),
        name = 'Layer')
plotReducedDim(
    sce_layer_sfac,
    dimred = 'PCA',
    colour_by = 'sample_name',
    theme_size = 20
)
plotReducedDim(
    sce_layer_sfac,
    dimred = 'PCA',
    colour_by = 'subject_position',
    theme_size = 20
)
plotReducedDim(
    sce_layer_sfac,
    dimred = 'PCA',
    colour_by = 'subject',
    theme_size = 20
)
plotReducedDim(
    sce_layer_sfac,
    dimred = 'PCA',
    colour_by = 'position',
    theme_size = 20
)
dev.off()
## Comparing the two, it seems to me that the default library size reduces the
## layer differences but makes more sense in the PCA plot:
## like it's PC1 for WM vs non-WM, then PC2 across layers 1 through 6
## instead of PC1 for some layers vs others, PC2 for WM vs layer 2?
## with layers like 6 and 5 across the PC1 vs PC2 plot.




pdf('pdf/reduced_dim_TSNE_perplexity5.pdf', useDingbats = FALSE)
plotReducedDim(
    sce_layer,
    dimred = 'TSNE_perplexity5',
    colour_by = 'layer_guess_reordered',
    theme_size = 20,
    point_size = 5
) +
    ggplot2::scale_fill_manual(values =  unname(Polychrome::palette36.colors(7))[c(2:7, 1)],
        name = 'Layer')
plotReducedDim(
    sce_layer,
    dimred = 'TSNE_perplexity5',
    colour_by = 'sample_name',
    theme_size = 20,
    point_size = 5
)
plotReducedDim(
    sce_layer,
    dimred = 'TSNE_perplexity5',
    colour_by = 'subject_position',
    theme_size = 20,
    point_size = 5
)
plotReducedDim(
    sce_layer,
    dimred = 'TSNE_perplexity5',
    colour_by = 'subject',
    theme_size = 20,
    point_size = 5
)
plotReducedDim(
    sce_layer,
    dimred = 'TSNE_perplexity5',
    colour_by = 'position',
    theme_size = 20,
    point_size = 5
)
dev.off()

pdf('pdf/reduced_dim_TSNE_perplexity15.pdf', useDingbats = FALSE)
plotReducedDim(
    sce_layer,
    dimred = 'TSNE_perplexity15',
    colour_by = 'layer_guess_reordered',
    theme_size = 20,
    point_size = 5
) +
    ggplot2::scale_fill_manual(values =  unname(Polychrome::palette36.colors(7))[c(2:7, 1)],
        name = 'Layer')
plotReducedDim(
    sce_layer,
    dimred = 'TSNE_perplexity15',
    colour_by = 'sample_name',
    theme_size = 20,
    point_size = 5
)
plotReducedDim(
    sce_layer,
    dimred = 'TSNE_perplexity15',
    colour_by = 'subject_position',
    theme_size = 20,
    point_size = 5
)
plotReducedDim(
    sce_layer,
    dimred = 'TSNE_perplexity15',
    colour_by = 'subject',
    theme_size = 20,
    point_size = 5
)
plotReducedDim(
    sce_layer,
    dimred = 'TSNE_perplexity15',
    colour_by = 'position',
    theme_size = 20,
    point_size = 5
)
dev.off()

pdf('pdf/reduced_dim_TSNE_perplexity20.pdf', useDingbats = FALSE)
plotReducedDim(
    sce_layer,
    dimred = 'TSNE_perplexity20',
    colour_by = 'layer_guess_reordered',
    theme_size = 20,
    point_size = 5
) +
    ggplot2::scale_fill_manual(values =  unname(Polychrome::palette36.colors(7))[c(2:7, 1)],
        name = 'Layer')
plotReducedDim(
    sce_layer,
    dimred = 'TSNE_perplexity20',
    colour_by = 'sample_name',
    theme_size = 20,
    point_size = 5
)
plotReducedDim(
    sce_layer,
    dimred = 'TSNE_perplexity20',
    colour_by = 'subject_position',
    theme_size = 20,
    point_size = 5
)
plotReducedDim(
    sce_layer,
    dimred = 'TSNE_perplexity20',
    colour_by = 'subject',
    theme_size = 20,
    point_size = 5
)
plotReducedDim(
    sce_layer,
    dimred = 'TSNE_perplexity20',
    colour_by = 'position',
    theme_size = 20,
    point_size = 5
)
dev.off()

pdf('pdf/reduced_dim_UMAP_neighbors15.pdf', useDingbats = FALSE)
plotReducedDim(
    sce_layer,
    dimred = 'UMAP_neighbors15',
    colour_by = 'layer_guess_reordered',
    theme_size = 20,
    point_size = 5
) +
    ggplot2::scale_fill_manual(values =  unname(Polychrome::palette36.colors(7))[c(2:7, 1)],
        name = 'Layer')
plotReducedDim(
    sce_layer,
    dimred = 'UMAP_neighbors15',
    colour_by = 'sample_name',
    theme_size = 20,
    point_size = 5
)
plotReducedDim(
    sce_layer,
    dimred = 'UMAP_neighbors15',
    colour_by = 'subject_position',
    theme_size = 20,
    point_size = 5
)
plotReducedDim(
    sce_layer,
    dimred = 'UMAP_neighbors15',
    colour_by = 'subject',
    theme_size = 20,
    point_size = 5
)
plotReducedDim(
    sce_layer,
    dimred = 'UMAP_neighbors15',
    colour_by = 'position',
    theme_size = 20,
    point_size = 5
)
dev.off()

## Across PCA, TSNE and UMAP the WM for samples 151669 and 151670
## (middle subject, position 0) group closer to the layer 6 data from
## the other images.
## In TSNE and UMAP, there points do group by subject more than they
## do in PCA space.


set.seed(20200121)
sced <-
    denoisePCA(
        sce_layer,
        dec,
        subset.row = top.hvgs,
        max.rank = 20,
        min.rank = 2
    )
ncol(reducedDim(sced, "PCA"))
# [1] 7

## Unlike here('Analysis', 'sce_scran.R') here this code ran within 2-5 secs
choices <-
    getClusteredPCs(reducedDim(sce_layer),
        max.rank = 20,
        min.rank = 2)
npcs <- metadata(choices)$chosen
npcs
# [1] 5
reducedDim(sce_layer, "PCAsub") <-
    reducedDim(sce_layer, "PCA")[, seq_len(npcs), drop = FALSE]

## This plot doesn't seem very informative in this scenario...
pdf('pdf/PC_choices.pdf', useDingbats = FALSE)
plot(choices$n.pcs,
    choices$n.clusters,
    xlab = "Number of PCs",
    ylab = "Number of clusters")
abline(a = 1, b = 1, col = "red")
abline(v = metadata(choices)$chosen,
    col = "grey80",
    lty = 2)
dev.off()

g_k5 <- buildSNNGraph(sce_layer, k = 5, use.dimred = 'PCA')
g_walk_k5 <- igraph::cluster_walktrap(g_k5)
clust_k5 <- sort_clusters(g_walk_k5$membership)
save(g_k5, g_walk_k5, file = 'rda/g_k5.Rdata')


addmargins(table(clust_k5, sce_layer$layer_guess))
# clust_k5 WM Layer 1 Layer 2 Layer 3 Layer 4 Layer 5 Layer 6 Sum
#      1    0       8       8      12       0       0       0  28
#      2    0       0       0       0      12      12       0  24
#      3    3       0       0       0       0       0      12  15
#      4    9       0       0       0       0       0       0   9
#      Sum 12       8       8      12      12      12      12  76

## ehem... that's why it makes more sense to have k = 7 :P
## since we have 7 layers. That was the k that denoisePCA()
## suggested earlier

g_k7 <- buildSNNGraph(sce_layer, k = 7, use.dimred = 'PCA')
g_walk_k7 <- igraph::cluster_walktrap(g_k7)
clust_k7 <- sort_clusters(g_walk_k7$membership)
save(g_k7, g_walk_k7, file = 'rda/g_k7.Rdata')

addmargins(table(clust_k7, sce_layer$layer_guess))
# clust_k7 WM Layer 1 Layer 2 Layer 3 Layer 4 Layer 5 Layer 6 Sum
#      1    0       8       8      12      10       0       0  38
#      2    4       0       0       0       2      12      12  30
#      3    8       0       0       0       0       0       0   8
#      Sum 12       8       8      12      12      12      12  76


## Hm... cut at k = 7
clust_k5_k7 <- sort_clusters(igraph::cut_at(g_walk_k5, n = 7))
clust_k7_k7 <- sort_clusters(igraph::cut_at(g_walk_k7, n = 7))

addmargins(table(clust_k5_k7, sce_layer$layer_guess))
# clust_k5_k7 WM Layer 1 Layer 2 Layer 3 Layer 4 Layer 5 Layer 6 Sum
#         1    0       0       8      12       0       0       0  20
#         2    0       0       0       0       4      12       0  16
#         3    3       0       0       0       0       0      12  15
#         4    9       0       0       0       0       0       0   9
#         5    0       8       0       0       0       0       0   8
#         6    0       0       0       0       4       0       0   4
#         7    0       0       0       0       4       0       0   4
#         Sum 12       8       8      12      12      12      12  76
addmargins(table(clust_k7_k7, sce_layer$layer_guess))
# clust_k7_k7 WM Layer 1 Layer 2 Layer 3 Layer 4 Layer 5 Layer 6 Sum
#         1    0       0       4      12       0       0       0  16
#         2    0       0       0       0       2      12       0  14
#         3    0       1       0       0      10       0       0  11
#         4    0       7       4       0       0       0       0  11
#         5    2       0       0       0       0       0       8  10
#         6    8       0       0       0       0       0       0   8
#         7    2       0       0       0       0       0       4   6
#         Sum 12       8       8      12      12      12      12  76




## Try with all 20 PCs
g_k20 <- buildSNNGraph(sce_layer, k = 20, use.dimred = 'PCA')
g_walk_k20 <- igraph::cluster_walktrap(g_k20)
clust_k20 <- sort_clusters(g_walk_k20$membership)
save(g_k20, g_walk_k20, file = 'rda/g_k20.Rdata')
clust_k20_k7 <- sort_clusters(igraph::cut_at(g_walk_k20, n = 7))

addmargins(table(clust_k20, sce_layer$layer_guess))
# clust_k20 WM Layer 1 Layer 2 Layer 3 Layer 4 Layer 5 Layer 6 Sum
#       1    0       8       8      12      12       0       0  40
#       2   12       0       0       0       0      12      12  36
#       Sum 12       8       8      12      12      12      12  76
addmargins(table(clust_k20_k7, sce_layer$layer_guess))
# clust_k20_k7 WM Layer 1 Layer 2 Layer 3 Layer 4 Layer 5 Layer 6 Sum
#          1    6       0       0       0       0       0      11  17
#          2    0       5       1      11       0       0       0  17
#          3    1       0       0       0       0      12       1  14
#          4    0       1       0       0      12       0       0  13
#          5    0       0       7       0       0       0       0   7
#          6    5       0       0       0       0       0       0   5
#          7    0       2       0       1       0       0       0   3
#          Sum 12       8       8      12      12      12      12  76

## They are all highly variable
addmargins(table(clust_k7_k7, clust_k5_k7))
#            clust_k5_k7
# clust_k7_k7  1  2  3  4  5  6  7 Sum
#         1   16  0  0  0  0  0  0  16
#         2    0 14  0  0  0  0  0  14
#         3    0  2  0  0  1  4  4  11
#         4    4  0  0  0  7  0  0  11
#         5    0  0 10  0  0  0  0  10
#         6    0  0  0  8  0  0  0   8
#         7    0  0  5  1  0  0  0   6
#         Sum 20 16 15  9  8  4  4  76

addmargins(table(clust_k7_k7, clust_k20_k7))
#            clust_k20_k7
# clust_k7_k7  1  2  3  4  5  6  7 Sum
#         1    0 12  0  0  3  0  1  16
#         2    0  0 12  2  0  0  0  14
#         3    0  0  0 11  0  0  0  11
#         4    0  5  0  0  4  0  2  11
#         5   10  0  0  0  0  0  0  10
#         6    3  0  0  0  0  5  0   8
#         7    4  0  2  0  0  0  0   6
#         Sum 17 17 14 13  7  5  3  76



## Hm... could it be that it's picking up the subject variability?
addmargins(table(clust_k7_k7, sce_layer$subject))
# clust_k7_k7 Br5292 Br5595 Br8100 Sum
#         1        4      4      8  16
#         2        4      6      4  14
#         3        4      2      5  11
#         4        8      0      3  11
#         5        4      6      0  10
#         6        4      0      4   8
#         7        0      2      4   6
#         Sum     28     20     28  76
addmargins(table(clust_k7_k7, sce_layer$layer_guess, sce_layer$subject))
# , ,  = Br5292
#
#
# clust_k7_k7 WM Layer 1 Layer 2 Layer 3 Layer 4 Layer 5 Layer 6 Sum
#         1    0       0       0       4       0       0       0   4
#         2    0       0       0       0       0       4       0   4
#         3    0       0       0       0       4       0       0   4
#         4    0       4       4       0       0       0       0   8
#         5    0       0       0       0       0       0       4   4
#         6    4       0       0       0       0       0       0   4
#         7    0       0       0       0       0       0       0   0
#         Sum  4       4       4       4       4       4       4  28
#
# , ,  = Br5595
#
#
# clust_k7_k7 WM Layer 1 Layer 2 Layer 3 Layer 4 Layer 5 Layer 6 Sum
#         1    0       0       0       4       0       0       0   4
#         2    0       0       0       0       2       4       0   6
#         3    0       0       0       0       2       0       0   2
#         4    0       0       0       0       0       0       0   0
#         5    2       0       0       0       0       0       4   6
#         6    0       0       0       0       0       0       0   0
#         7    2       0       0       0       0       0       0   2
#         Sum  4       0       0       4       4       4       4  20
#
# , ,  = Br8100
#
#
# clust_k7_k7 WM Layer 1 Layer 2 Layer 3 Layer 4 Layer 5 Layer 6 Sum
#         1    0       0       4       4       0       0       0   8
#         2    0       0       0       0       0       4       0   4
#         3    0       1       0       0       4       0       0   5
#         4    0       3       0       0       0       0       0   3
#         5    0       0       0       0       0       0       0   0
#         6    4       0       0       0       0       0       0   4
#         7    0       0       0       0       0       0       4   4
#         Sum  4       4       4       4       4       4       4  28
#
# , ,  = Sum
#
#
# clust_k7_k7 WM Layer 1 Layer 2 Layer 3 Layer 4 Layer 5 Layer 6 Sum
#         1    0       0       4      12       0       0       0  16
#         2    0       0       0       0       2      12       0  14
#         3    0       1       0       0      10       0       0  11
#         4    0       7       4       0       0       0       0  11
#         5    2       0       0       0       0       0       8  10
#         6    8       0       0       0       0       0       0   8
#         7    2       0       0       0       0       0       4   6
#         Sum 12       8       8      12      12      12      12  76

## So subject 1, it merges layers 1 & 2
## subject 2, it breaks layer 4 (mixes it a bit with layer 5) and WM
## subject 3, for one it merges layer 1 with 4, merges layers 2 and 3



addmargins(table(clust_k5_k7, sce_layer$subject))
# clust_k5_k7 Br5292 Br5595 Br8100 Sum
#         1        8      4      8  20
#         2        4      8      4  16
#         3        4      7      4  15
#         4        4      1      4   9
#         5        4      0      4   8
#         6        4      0      0   4
#         7        0      0      4   4
#         Sum     28     20     28  76

addmargins(table(clust_k5_k7, sce_layer$layer_guess, sce_layer$subject))
# , ,  = Br5292
#
#
# clust_k5_k7 WM Layer 1 Layer 2 Layer 3 Layer 4 Layer 5 Layer 6 Sum
#         1    0       0       4       4       0       0       0   8
#         2    0       0       0       0       0       4       0   4
#         3    0       0       0       0       0       0       4   4
#         4    4       0       0       0       0       0       0   4
#         5    0       4       0       0       0       0       0   4
#         6    0       0       0       0       4       0       0   4
#         7    0       0       0       0       0       0       0   0
#         Sum  4       4       4       4       4       4       4  28
#
# , ,  = Br5595
#
#
# clust_k5_k7 WM Layer 1 Layer 2 Layer 3 Layer 4 Layer 5 Layer 6 Sum
#         1    0       0       0       4       0       0       0   4
#         2    0       0       0       0       4       4       0   8
#         3    3       0       0       0       0       0       4   7
#         4    1       0       0       0       0       0       0   1
#         5    0       0       0       0       0       0       0   0
#         6    0       0       0       0       0       0       0   0
#         7    0       0       0       0       0       0       0   0
#         Sum  4       0       0       4       4       4       4  20
#
# , ,  = Br8100
#
#
# clust_k5_k7 WM Layer 1 Layer 2 Layer 3 Layer 4 Layer 5 Layer 6 Sum
#         1    0       0       4       4       0       0       0   8
#         2    0       0       0       0       0       4       0   4
#         3    0       0       0       0       0       0       4   4
#         4    4       0       0       0       0       0       0   4
#         5    0       4       0       0       0       0       0   4
#         6    0       0       0       0       0       0       0   0
#         7    0       0       0       0       4       0       0   4
#         Sum  4       4       4       4       4       4       4  28
#
# , ,  = Sum
#
#
# clust_k5_k7 WM Layer 1 Layer 2 Layer 3 Layer 4 Layer 5 Layer 6 Sum
#         1    0       0       8      12       0       0       0  20
#         2    0       0       0       0       4      12       0  16
#         3    3       0       0       0       0       0      12  15
#         4    9       0       0       0       0       0       0   9
#         5    0       8       0       0       0       0       0   8
#         6    0       0       0       0       4       0       0   4
#         7    0       0       0       0       4       0       0   4
#         Sum 12       8       8      12      12      12      12  76

## In this one, again for subject 1, it merges layers 1 & 2
## subject 2, it breaks WM, merges layers 4 and 5, WM and layer 6
## subject 3, for one it merges layers 2 and 3


## Try some of the k-means clustering
set.seed(20200122)
clust.kmeans <- kmeans(reducedDim(sce_layer, "PCA"), centers = 7)
save(clust.kmeans, file = 'rda/clust.kmeans.Rdata')
table(sort_clusters(clust.kmeans$cluster))
#  1  2  3  4  5  6  7
# 22 16 14  9  7  4  4

addmargins(table(
    sort_clusters(clust.kmeans$cluster),
    sce_layer$layer_guess,
    sce_layer$subject
))
# , ,  = Br5292
#
#
#       WM Layer1 Layer2 Layer3 Layer4 Layer5 Layer6 Sum
#   1    0      0      0      4      4      0      0   8
#   2    0      4      4      0      0      0      0   8
#   3    0      0      0      0      0      4      0   4
#   4    0      0      0      0      0      0      3   3
#   5    3      0      0      0      0      0      0   3
#   6    0      0      0      0      0      0      0   0
#   7    1      0      0      0      0      0      1   2
#   Sum  4      4      4      4      4      4      4  28
#
# , ,  = Br5595
#
#
#       WM Layer1 Layer2 Layer3 Layer4 Layer5 Layer6 Sum
#   1    0      0      0      4      2      0      0   6
#   2    0      0      0      0      0      0      0   0
#   3    0      0      0      0      2      4      0   6
#   4    2      0      0      0      0      0      4   6
#   5    0      0      0      0      0      0      0   0
#   6    0      0      0      0      0      0      0   0
#   7    2      0      0      0      0      0      0   2
#   Sum  4      0      0      4      4      4      4  20
#
# , ,  = Br8100
#
#
#       WM Layer1 Layer2 Layer3 Layer4 Layer5 Layer6 Sum
#   1    0      0      0      4      4      0      0   8
#   2    0      4      4      0      0      0      0   8
#   3    0      0      0      0      0      4      0   4
#   4    0      0      0      0      0      0      0   0
#   5    4      0      0      0      0      0      0   4
#   6    0      0      0      0      0      0      4   4
#   7    0      0      0      0      0      0      0   0
#   Sum  4      4      4      4      4      4      4  28
#
# , ,  = Sum
#
#
#       WM Layer1 Layer2 Layer3 Layer4 Layer5 Layer6 Sum
#   1    0      0      0     12     10      0      0  22
#   2    0      8      8      0      0      0      0  16
#   3    0      0      0      0      2     12      0  14
#   4    2      0      0      0      0      0      7   9
#   5    7      0      0      0      0      0      0   7
#   6    0      0      0      0      0      0      4   4
#   7    3      0      0      0      0      0      1   4
#   Sum 12      8      8     12     12     12     12  76

set.seed(20200122)
gaps <- clusGap(reducedDim(sce_layer, "PCA"), kmeans, K.max = 20)
best.k <- maxSE(gaps$Tab[, "gap"], gaps$Tab[, "SE.sim"])
best.k
# [1] 7
## Yay, at least here we get 7 :P

ncells <- tabulate(sort_clusters(clust.kmeans$cluster))
tab <- data.frame(wcss = clust.kmeans$withinss, ncells = ncells)
tab$rms <- sqrt(tab$wcss / tab$ncells)
tab
#         wcss ncells       rms
# 1  2501.3521     22 10.662918
# 2  4071.9149     16 15.952890
# 3  5505.8851     14 19.831225
# 4   246.5383      9  5.233846
# 5  9272.5498      7 36.395741
# 6 10386.1113      4 50.956136
# 7  5770.5179      4 37.981963

pdf('pdf/kmeans_gaps.pdf', useDingbats = FALSE)
plot(gaps$Tab[, "gap"], xlab = "Number of clusters", ylab = "Gap statistic")
abline(v = best.k, col = "red")
dev.off()

centers <- clust.kmeans$centers
## Fix the names
rownames(centers) <- unique(sort_clusters(clust.kmeans$cluster))
cent.tree <- hclust(dist(centers), "ward.D2")
pdf('pdf/kmeans_k7_dendro.pdf', useDingbats = FALSE)
plot(cent.tree)
dev.off()



## Save for later
colData(sce_layer)$c_k5_k7 <- clust_k5_k7
colData(sce_layer)$c_k7_k7 <- clust_k7_k7
colData(sce_layer)$c_k20_k7 <- clust_k20_k7
colData(sce_layer)$kmeans_k7 <- sort_clusters(clust.kmeans$cluster)

save(sce_layer, top.hvgs, file = 'rda/sce_layer.Rdata')
## For mapping back to the original sce object
save(layerIndexes, file = 'rda/layerIndexes.Rdata')
## For subsetting again the genes if necessary
save(ix_mito, selected_genes, file = 'rda/selected_genes.Rdata')




## Try using the scran::findMarkers() function instead of coding the
## looping code myself around limma::lmFit()
## More at https://osca.bioconductor.org/marker-detection.html

## Use 'block' for random effects
## https://support.bioconductor.org/p/29768/

## I see that there's more info about findMarkers() now at
## https://bioconductor.riken.jp/packages/release/workflows/vignettes/simpleSingleCell/inst/doc/de.html
## and
## https://bioconductor.org/packages/3.9/workflows/vignettes/simpleSingleCell/inst/doc/reads.html#detecting-marker-genes-between-clusters
markers_layer <- lapply(c('any', 'all'), function(pval) {
    directions <- c('any', 'up')
    res_direc <- lapply(directions, function(direc) {
        message(paste(
            Sys.time(),
            'finding markers for p-val',
            pval,
            'and',
            direc,
            'direction'
        ))
        findMarkers(
            sce_layer,
            sce_layer$layer_guess,
            pval.type = pval,
            direction = direc,
            block = sce_layer$subject_position,
            gene.names = rowData(sce_layer)$gene_name,
            full.stats = FALSE ## Setting this to TRUE messes up the plotting code right now
        )
    })
    names(res_direc) <- directions
    return(res_direc)
})
names(markers_layer) <- c('any', 'all')
save(markers_layer, file = 'rda/markers_layer.Rdata')



markers_layer_wilcox <- lapply(c('any', 'all'), function(pval) {
    directions <- c('any', 'up')
    res_direc <- lapply(directions, function(direc) {
        message(paste(
            Sys.time(),
            'finding markers for p-val',
            pval,
            'and',
            direc,
            'direction'
        ))
        findMarkers(
            sce_layer,
            sce_layer$layer_guess,
            pval.type = pval,
            direction = direc,
            block = sce_layer$subject_position,
            gene.names = rowData(sce_layer)$gene_name,
            test.type = 'wilcox',
            full.stats = FALSE
        )
    })
    names(res_direc) <- directions
    return(res_direc)
})
names(markers_layer_wilcox) <- c('any', 'all')
save(markers_layer_wilcox, file = 'rda/markers_layer_wilcox.Rdata')





summary(gene_ann(rowData(sce_layer)$gene_name))
#  KM_Zeng          BM         RNAscope
# FALSE:22255   FALSE:22275   FALSE:18594
# TRUE :   76   TRUE :   56   TRUE : 3737
nrow(genes_km_raw)
# [1] 81
nrow(genes_bm_raw)
# [1] 65
## So 5 and 9 genes are not here to begin with

summary(gene_ann(rowData(sce_original)$gene_name))
#  KM_Zeng          BM         RNAscope
# FALSE:33461   FALSE:33481   FALSE:29453
# TRUE :   77   TRUE :   57   TRUE : 4085
## But only 1 of each wasn't in the original data

## Find which are these genes:
subset(gene_ann(rowData(sce_original)$gene_name[-ix_mito][-selected_genes]), KM_Zeng == 'TRUE' |
        BM == 'TRUE')
#        KM_Zeng    BM
# LHX5     FALSE  TRUE
# CACNG5    TRUE FALSE

with(gene_ann(rowData(sce_original)$gene_name), addmargins(table(KM_Zeng, BM)))
#        BM
# KM_Zeng FALSE  TRUE   Sum
#   FALSE 33416    45 33461
#   TRUE     65    12    77
#   Sum   33481    57 33538

## Find which ones are absent to begin with
genes_km_raw$Gene[!tolower(genes_km_raw$Gene) %in% tolower(rowData(sce_original)$gene_name)]
# [1] "C4orf31"   "Pvrl3"     "C20orf103" "Kiaa1456"
genes_bm_raw$Gene[!tolower(genes_bm_raw$Gene) %in% tolower(rowData(sce_original)$gene_name)]
# [1] "A930038C07Rik" "9830123M21Rik" "Brn1"          "Brn2"
# [5] "6430573F11Rik" "C030003D03Rik" "Trb"           "Igh6"

## Double check
rowData(sce_original)$gene_name[grep('c4orf3', tolower(rowData(sce_original)$gene_name))]
rowData(sce_original)$gene_name[grep('c20orf', tolower(rowData(sce_original)$gene_name))]
rowData(sce_original)$gene_name[grep('pvr', tolower(rowData(sce_original)$gene_name))]
rowData(sce_original)$gene_name[grep('rik', tolower(rowData(sce_original)$gene_name))]
# [1] "GRIK3"     "GRIK2"     "GRIK4"     "GRIK5"     "GRIK1"     "GRIK1-AS1"




## t-test versions
plot_markers_loop(
    markers_layer,
    'markers_t-test_logFC',
    plot_markers_logfc,
    prefix = 'logFC',
    breaks = seq(-2, 2, length.out = 101)
)
plot_markers_loop(
    markers_layer,
    'markers_t-test_expr',
    plot_markers_expr,
    h = 17,
    color = viridis::viridis(101)
)
plot_markers_loop(
    markers_layer,
    'markers_t-test_expr_centered',
    plot_markers_expr,
    h = 17,
    center = TRUE,
    zlim = c(-2, 2)
)

## wilcox test versions
plot_markers_loop(
    markers_layer_wilcox,
    'markers_wilcox_AUC',
    plot_markers_logfc,
    prefix = 'AUC',
    breaks = seq(0, 1, length.out = 101),
    color = viridis::viridis(101)
)
plot_markers_loop(
    markers_layer_wilcox,
    'markers_wilcox_expr',
    plot_markers_expr,
    h = 17,
    color = viridis::viridis(101)
)
plot_markers_loop(
    markers_layer_wilcox,
    'markers_wilcox_expr_centered',
    plot_markers_expr,
    h = 17,
    center = TRUE,
    zlim = c(-2, 2)
)



find_marker_gene('MOBP', markers_layer, 'any', 'any')
find_marker_gene('MOBP', markers_layer, 'any', 'up')
find_marker_gene('MOBP', markers_layer, 'all', 'any')
find_marker_gene('MOBP', markers_layer, 'all', 'up')







###### Direct limma approach ####
#################################


## Extract the data
mat <- assays(sce_layer)$logcounts

## Build a group model
mod <- with(colData(sce_layer), model.matrix(~ 0 + layer_guess))
colnames(mod) <- gsub('layer_guess', '', colnames(mod))
## Takes like 2 min to run
corfit <-
    duplicateCorrelation(mat, mod, block = sce_layer$subject_position)
fit <-
    lmFit(
        mat,
        design = mod,
        block = sce_layer$subject_position,
        correlation = corfit$consensus.correlation
    )
eb <- eBayes(fit)


## Define the contrasts for each layer vs the rest (excluding WM comparisons since we have that one already)
layer_combs <- combn(colnames(mod), 2)
layer_contrasts <- apply(layer_combs, 2, function(x) {
    z <- paste(x, collapse = '-')
    makeContrasts(contrasts = z, levels = mod)
})
rownames(layer_contrasts) <- colnames(mod)
colnames(layer_contrasts) <-
    apply(layer_combs, 2, paste, collapse = '-')
eb_contrasts <- eBayes(contrasts.fit(fit, layer_contrasts))


## Extract the p-values and add the WM comparisons too
pvals_contrasts <- eb_contrasts$p.value

## Fing the sig ones
data.frame(
    'FDRsig' = colSums(apply(pvals_contrasts, 2, p.adjust, 'fdr') < 0.05),
    'Pval10-6sig' = colSums(pvals_contrasts < 1e-6),
    'Pval10-8sig' = colSums(pvals_contrasts < 1e-8)
)
#               FDRsig Pval10.6sig Pval10.8sig
# WM-Layer1       5667        1291         766
# WM-Layer2       7992        2529        1652
# WM-Layer3       8500        2875        1986
# WM-Layer4       8458        2874        1944
# WM-Layer5       7979        2653        1752
# WM-Layer6       6686        1922        1210
# Layer1-Layer2   3686         486         185
# Layer1-Layer3   3566         641         267
# Layer1-Layer4   4677         940         512
# Layer1-Layer5   4638        1090         599
# Layer1-Layer6   4255         893         459
# Layer2-Layer3    377          70          24
# Layer2-Layer4   2284         415         213
# Layer2-Layer5   2305         508         275
# Layer2-Layer6   2447         467         256
# Layer3-Layer4    327          60          24
# Layer3-Layer5    942         224         124
# Layer3-Layer6   1752         334         192
# Layer4-Layer5    292          67          35
# Layer4-Layer6   1745         246         140
# Layer5-Layer6    515         103          44






## Next for each layer test that layer vs the rest
layer_idx <- splitit(sce_layer$layer_guess)

eb0_list <- lapply(layer_idx, function(x) {
    res <- rep(0, ncol(sce_layer))
    res[x] <- 1
    m <- model.matrix(~ res)
    eBayes(
        lmFit(
            mat,
            design = m,
            block = sce_layer$subject_position,
            correlation = corfit$consensus.correlation
        )
    )
})

## Extract the p-values
pvals0_contrasts <- sapply(eb0_list, function(x) {
    x$p.value[, 2, drop = FALSE]
})

data.frame(
    'FDRsig' = colSums(apply(pvals0_contrasts, 2, p.adjust, 'fdr') < 0.05),
    'Pval10-6sig' = colSums(pvals0_contrasts < 1e-6),
    'Pval10-8sig' = colSums(pvals0_contrasts < 1e-8)
)
#        FDRsig Pval10.6sig Pval10.8sig
# WM       9124        3039        2021
# Layer1   3033         368         170
# Layer2   1562         128          28
# Layer3    183           9           3
# Layer4    740          29          10
# Layer5    643          64          27
# Layer6    379          71          27

## Extract the tstats
tstats0_contrasts <- sapply(eb0_list, function(x) {
    x$t[, 2, drop = FALSE]
})

## Save for later
save(eb0_list, file = 'rda/eb0_list.Rdata')
save(eb_contrasts, file = 'rda/eb_contrasts.Rdata')


## Extract data
sig_genes_layer0 <-
    sig_genes_extract(tstats0_contrasts, pvals0_contrasts)

sig_genes_summary <- function(x) {
    summary(x[!duplicated(x$ensembl), c('KM_Zeng', 'BM', 'RNAscope')])
}

sig_genes_summary(sig_genes_layer0)
#  KM_Zeng       BM      RNAscope
# FALSE:62   FALSE:66   FALSE:35
# TRUE : 8   TRUE : 4   TRUE :35



## Now extract the data for the paired comparisons
tstats_contrasts <- eb_contrasts$t
sig_genes_layer <-
    sig_genes_extract(tstats_contrasts, pvals_contrasts)
sig_genes_summary(sig_genes_layer)
#  KM_Zeng       BM      RNAscope
# FALSE:80   FALSE:84   FALSE:44
# TRUE :11   TRUE : 7   TRUE :47


## likely best to extract
sig_genes_layer_rev <-
    sig_genes_extract(-1 * tstats_contrasts, pvals_contrasts)
sig_genes_layer_rev$layer <-
    rep(apply(layer_combs[c(2, 1),], 2, paste, collapse = '-'), each = 10)
sig_genes_summary(sig_genes_layer_rev)
# KM_Zeng       BM      RNAscope
# FALSE:78   FALSE:90   FALSE:51
# TRUE :24   TRUE :12   TRUE :51

sig_genes <- rbind(
    cbind(sig_genes_layer0, test = 'layer_vs_rest'),
    cbind(sig_genes_layer, test = 'paired_layers'),
    cbind(sig_genes_layer_rev, test = 'paired_layers')
)
sig_genes_summary(sig_genes)
#  KM_Zeng        BM       RNAscope
# FALSE:170   FALSE:185   FALSE:100
# TRUE : 28   TRUE : 13   TRUE : 98

sig_genes_unique <- splitit(sig_genes$ensembl)
length(sig_genes_unique)
# [1] 198

sig_genes <- DataFrame(sig_genes)
sig_genes$in_rows <-
    IntegerList(sig_genes_unique)[sig_genes$ensembl]
sig_genes$results <-
    CharacterList(sapply(sig_genes$in_rows, function(x)
        paste0(sig_genes$layer[x], '_top', sig_genes$top[x])))


sig_genes_df <- sig_genes
## Fix for writing to csv
sig_genes_df$in_rows <-
    sapply(sig_genes_df$in_rows, paste0, collapse = ';')
sig_genes_df$results <-
    sapply(sig_genes_df$results, paste0, collapse = ';')
write.csv(sig_genes_df,
    file = 'sig_genes.csv',
    quote = FALSE,
    row.names = FALSE)

## Save lists
save(sig_genes,
    sig_genes_layer0,
    sig_genes_layer,
    sig_genes_layer_rev,
    file = 'rda/layer_sig_genes.Rdata')


## Make gene grid plots
pdf_dir <- 'pdf/gene_grid/sig_genes'
dir.create(pdf_dir, showWarnings = FALSE, recursive = TRUE)


## Only make the plots for the unique ones
## Takes about 1 - 1.5 hours
for (i in match(names(sig_genes_unique), sig_genes$ensembl)) {
    # i <- 1
    message(paste(Sys.time(), 'making the plot for', i, 'gene', sig_genes$gene[i]))
    sce_image_grid_gene(
        sce_original,
        geneid = paste0(sig_genes$gene[i], '; ', sig_genes$ensembl[i]),
        pdf_file = file.path(pdf_dir,
            paste0(
                sig_genes$gene[i],
                '_',
                gsub('top', 'r', gsub('Layer', 'L', sig_genes_df$results[i])),
                '.pdf'
            )),
        ... = gsub('top', 'r', gsub('Layer', 'L', sig_genes_df$results[i]))
    )
}

## Make heatmaps
rownames(sig_genes_df) <- sig_genes_df$gene

pdf('pdf/markers_layer_expr.pdf',
    useDingbats = FALSE,
    height = 17)
xx <- plot_markers_expr(
    split(sig_genes_df, sig_genes$layer),
    pval.type = 'all',
    color = viridis::viridis(101)
)
dev.off()

pdf('pdf/markers_layer_expr_centered.pdf',
    useDingbats = FALSE,
    height = 17)
xx <- plot_markers_expr(
    split(sig_genes_df, sig_genes$layer),
    pval.type = 'all',
    center = TRUE,
    zlim = c(-2, 2)
)
dev.off()

## Make boxplots
layer_guess_reordered <-
    factor(sce_layer$layer_guess, levels = c(paste0('Layer', 1:6), 'WM'))
pdf('pdf/markers_layer_boxplots.pdf', useDingbats = FALSE)
set.seed(20200206)
for (i in seq_len(nrow(sig_genes))) {
    # i <- 1
    message(paste(Sys.time(), 'making the plot for', i, 'gene', sig_genes$gene[i]))
    boxplot(
        mat[sig_genes$gene_index[i], ] ~ layer_guess_reordered,
        xlab = 'Layer',
        ylab = 'logcounts',
        main = paste(
            sig_genes$gene[i],
            sig_genes$ensembl[i],
            sig_genes$layer[i],
            '\n',
            
            'tstat',
            formatC(sig_genes$tstat[i], format = "e", digits = 2),
            'p',
            formatC(sig_genes$pval[i], format = "e", digits = 2),
            
            '\n',
            gsub('top', 'r', gsub('Layer', 'L', sig_genes_df$results[i]))
        ),
        outline = FALSE,
        cex = 1.5
    )
    points(
        mat[sig_genes$gene_index[i], ] ~ jitter(as.integer(layer_guess_reordered)),
        pch = 21,
        bg = Polychrome::palette36.colors(7)[as.integer(sce_layer$layer_guess)],
        cex = 1.5
    )
    # legend(
    #     "top",
    #     levels(sce_layer$layer_guess),
    #     col =  Polychrome::palette36.colors(7),
    #     lwd = 2
    # )
}
dev.off()


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────
#  setting  value
#  version  R version 3.6.1 Patched (2019-10-31 r77350)
#  os       CentOS Linux 7 (Core)
#  system   x86_64, linux-gnu
#  ui       X11
#  language (EN)
#  collate  en_US.UTF-8
#  ctype    en_US.UTF-8
#  tz       US/Eastern
#  date     2020-01-24
#
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────
#  package              * version   date       lib source
#  assertthat             0.2.1     2019-03-21 [2] CRAN (R 3.6.1)
#  backports              1.1.5     2019-10-02 [1] CRAN (R 3.6.1)
#  beeswarm               0.2.3     2016-04-25 [2] CRAN (R 3.6.1)
#  Biobase              * 2.46.0    2019-10-29 [2] Bioconductor
#  BiocGenerics         * 0.32.0    2019-10-29 [1] Bioconductor
#  BiocNeighbors          1.4.1     2019-11-01 [2] Bioconductor
#  BiocParallel         * 1.20.1    2019-12-21 [1] Bioconductor
#  BiocSingular           1.2.0     2019-10-29 [2] Bioconductor
#  bitops                 1.0-6     2013-08-17 [2] CRAN (R 3.6.1)
#  cellranger             1.1.0     2016-07-27 [1] CRAN (R 3.6.1)
#  cli                    2.0.0     2019-12-09 [1] CRAN (R 3.6.1)
#  cluster              * 2.1.0     2019-06-19 [3] CRAN (R 3.6.1)
#  colorout             * 1.2-2     2019-10-31 [1] Github (jalvesaq/colorout@641ed38)
#  colorspace             1.4-1     2019-03-18 [2] CRAN (R 3.6.1)
#  cowplot              * 1.0.0     2019-07-11 [1] CRAN (R 3.6.1)
#  crayon                 1.3.4     2017-09-16 [1] CRAN (R 3.6.1)
#  DelayedArray         * 0.12.0    2019-10-29 [2] Bioconductor
#  DelayedMatrixStats     1.8.0     2019-10-29 [2] Bioconductor
#  digest                 0.6.23    2019-11-23 [1] CRAN (R 3.6.1)
#  dplyr                  0.8.3     2019-07-04 [1] CRAN (R 3.6.1)
#  dqrng                  0.2.1     2019-05-17 [1] CRAN (R 3.6.1)
#  edgeR                  3.28.0    2019-10-29 [1] Bioconductor
#  fansi                  0.4.0     2018-10-05 [1] CRAN (R 3.6.1)
#  GenomeInfoDb         * 1.22.0    2019-10-29 [1] Bioconductor
#  GenomeInfoDbData       1.2.2     2019-10-28 [2] Bioconductor
#  GenomicRanges        * 1.38.0    2019-10-29 [1] Bioconductor
#  ggbeeswarm             0.6.0     2017-08-07 [2] CRAN (R 3.6.1)
#  ggplot2              * 3.2.1     2019-08-10 [1] CRAN (R 3.6.1)
#  glue                   1.3.1     2019-03-12 [1] CRAN (R 3.6.1)
#  googledrive            1.0.0     2019-08-19 [1] CRAN (R 3.6.1)
#  gridExtra              2.3       2017-09-09 [2] CRAN (R 3.6.1)
#  gtable                 0.3.0     2019-03-25 [2] CRAN (R 3.6.1)
#  here                 * 0.1       2017-05-28 [1] CRAN (R 3.6.1)
#  htmltools              0.4.0     2019-10-04 [1] CRAN (R 3.6.1)
#  htmlwidgets            1.5.1     2019-10-08 [1] CRAN (R 3.6.1)
#  httpuv                 1.5.2     2019-09-11 [1] CRAN (R 3.6.1)
#  igraph                 1.2.4.1   2019-04-22 [2] CRAN (R 3.6.1)
#  IRanges              * 2.20.1    2019-11-20 [1] Bioconductor
#  irlba                  2.3.3     2019-02-05 [1] CRAN (R 3.6.1)
#  jaffelab             * 0.99.29   2019-11-04 [1] Github (LieberInstitute/jaffelab@a7d87cb)
#  jsonlite               1.6       2018-12-07 [2] CRAN (R 3.6.1)
#  later                  1.0.0     2019-10-04 [1] CRAN (R 3.6.1)
#  lattice                0.20-38   2018-11-04 [3] CRAN (R 3.6.1)
#  lazyeval               0.2.2     2019-03-15 [2] CRAN (R 3.6.1)
#  limma                * 3.42.0    2019-10-29 [1] Bioconductor
#  locfit                 1.5-9.1   2013-04-20 [2] CRAN (R 3.6.1)
#  magrittr               1.5       2014-11-22 [1] CRAN (R 3.6.1)
#  Matrix                 1.2-17    2019-03-22 [3] CRAN (R 3.6.1)
#  matrixStats          * 0.55.0    2019-09-07 [1] CRAN (R 3.6.1)
#  munsell                0.5.0     2018-06-12 [2] CRAN (R 3.6.1)
#  pheatmap             * 1.0.12    2019-01-04 [2] CRAN (R 3.6.1)
#  pillar                 1.4.3     2019-12-20 [1] CRAN (R 3.6.1)
#  pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 3.6.1)
#  png                    0.1-7     2013-12-03 [2] CRAN (R 3.6.1)
#  Polychrome           * 1.2.3     2019-08-01 [1] CRAN (R 3.6.1)
#  promises               1.1.0     2019-10-04 [1] CRAN (R 3.6.1)
#  purrr                  0.3.3     2019-10-18 [2] CRAN (R 3.6.1)
#  R6                     2.4.0     2019-02-14 [2] CRAN (R 3.6.1)
#  rafalib              * 1.0.0     2015-08-09 [1] CRAN (R 3.6.1)
#  RColorBrewer           1.1-2     2014-12-07 [2] CRAN (R 3.6.1)
#  Rcpp                   1.0.3     2019-11-08 [1] CRAN (R 3.6.1)
#  RCurl                  1.95-4.12 2019-03-04 [2] CRAN (R 3.6.1)
#  readxl               * 1.3.1     2019-03-13 [2] CRAN (R 3.6.1)
#  rlang                  0.4.2     2019-11-23 [1] CRAN (R 3.6.1)
#  rmote                * 0.3.4     2019-10-31 [1] Github (cloudyr/rmote@fbce611)
#  rprojroot              1.3-2     2018-01-03 [2] CRAN (R 3.6.1)
#  rsvd                   1.0.2     2019-07-29 [1] CRAN (R 3.6.1)
#  S4Vectors            * 0.24.1    2019-12-01 [1] Bioconductor
#  scales                 1.0.0     2018-08-09 [2] CRAN (R 3.6.1)
#  scater               * 1.14.3    2019-11-07 [2] Bioconductor
#  scatterplot3d          0.3-41    2018-03-14 [1] CRAN (R 3.6.1)
#  scran                * 1.14.5    2019-11-19 [1] Bioconductor
#  segmented              1.0-0     2019-06-17 [2] CRAN (R 3.6.1)
#  servr                  0.15      2019-08-07 [1] CRAN (R 3.6.1)
#  sessioninfo          * 1.1.1     2018-11-05 [1] CRAN (R 3.6.1)
#  SingleCellExperiment * 1.8.0     2019-10-29 [2] Bioconductor
#  statmod                1.4.32    2019-05-29 [2] CRAN (R 3.6.1)
#  SummarizedExperiment * 1.16.1    2019-12-19 [1] Bioconductor
#  tibble                 2.1.3     2019-06-06 [1] CRAN (R 3.6.1)
#  tidyselect             0.2.5     2018-10-11 [2] CRAN (R 3.6.1)
#  vipor                  0.4.5     2017-03-22 [2] CRAN (R 3.6.1)
#  viridis                0.5.1     2018-03-29 [2] CRAN (R 3.6.1)
#  viridisLite          * 0.3.0     2018-02-01 [2] CRAN (R 3.6.1)
#  withr                  2.1.2     2018-03-15 [2] CRAN (R 3.6.1)
#  xfun                   0.11      2019-11-12 [1] CRAN (R 3.6.1)
#  XVector                0.26.0    2019-10-29 [1] Bioconductor
#  zlibbioc               1.32.0    2019-10-29 [2] Bioconductor
#
# [1] /users/lcollado/R/3.6.x
# [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/site-library
# [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/library

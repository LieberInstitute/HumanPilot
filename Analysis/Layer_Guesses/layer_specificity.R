library('SingleCellExperiment')
library('here')
library('jaffelab')
library('scater')
library('scran')
library('pheatmap')
library('readxl')
library('Polychrome')
library('sessioninfo')

dir.create('pdf', showWarnings = FALSE)
dir.create('rda', showWarnings = FALSE)

## Load data
load(here(
    'Analysis',
    'Human_DLPFC_Visium_processedData_sce_scran.Rdata'
))

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
sce$layer_guess <-
    factor(sce$layer_guess, levels = c('WM', paste('Layer', 1:6)))

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

## From https://gist.githubusercontent.com/mages/5339689/raw/2aaa482dfbbecbfcb726525a3d81661f9d802a8e/add.alpha.R
add.alpha <- function(col, alpha = 1) {
    if (missing(col))
        stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb) / 255, 2,
        function(x)
            rgb(x[1], x[2], x[3], alpha = alpha))
}

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


pdf('pdf/reduced_dim_PCA.pdf', useDingbats = FALSE)
plotReducedDim(sce_layer, dimred = 'PCA', colour_by = 'layer_guess') +
    ggplot2::scale_fill_manual(values =  unname(Polychrome::palette36.colors(7)),
        name = 'Layer')
plotReducedDim(sce_layer, dimred = 'PCA', colour_by = 'sample_name')
plotReducedDim(sce_layer, dimred = 'PCA', colour_by = 'subject_position')
plotReducedDim(sce_layer, dimred = 'PCA', colour_by = 'subject')
plotReducedDim(sce_layer, dimred = 'PCA', colour_by = 'position')
dev.off()

## Just to compare the PCA on the data using library size factors
## by sample repeated across layers
sce_layer_sfac <-
    runPCA(sce_layer_sfac,
        subset_row = top.hvgs,
        ncomponents = 20)
pdf('pdf/reduced_dim_PCA_sfact.pdf', useDingbats = FALSE)
plotReducedDim(sce_layer_sfac, dimred = 'PCA', colour_by = 'layer_guess') +
    ggplot2::scale_fill_manual(values =  unname(Polychrome::palette36.colors(7)),
        name = 'Layer')
plotReducedDim(sce_layer_sfac, dimred = 'PCA', colour_by = 'sample_name')
plotReducedDim(sce_layer_sfac, dimred = 'PCA', colour_by = 'subject_position')
plotReducedDim(sce_layer_sfac, dimred = 'PCA', colour_by = 'subject')
plotReducedDim(sce_layer_sfac, dimred = 'PCA', colour_by = 'position')
dev.off()
## Comparing the two, it seems to me that the default library size reduces the
## layer differences but makes more sense in the PCA plot:
## like it's PC1 for WM vs non-WM, then PC2 across layers 1 through 6
## instead of PC1 for some layers vs others, PC2 for WM vs layer 2?
## with layers like 6 and 5 across the PC1 vs PC2 plot.




pdf('pdf/reduced_dim_TSNE_perplexity5.pdf', useDingbats = FALSE)
plotReducedDim(sce_layer, dimred = 'TSNE_perplexity5', colour_by = 'layer_guess') +
    ggplot2::scale_fill_manual(values =  unname(Polychrome::palette36.colors(7)),
        name = 'Layer')
plotReducedDim(sce_layer, dimred = 'TSNE_perplexity5', colour_by = 'sample_name')
plotReducedDim(sce_layer, dimred = 'TSNE_perplexity5', colour_by = 'subject_position')
plotReducedDim(sce_layer, dimred = 'TSNE_perplexity5', colour_by = 'subject')
plotReducedDim(sce_layer, dimred = 'TSNE_perplexity5', colour_by = 'position')
dev.off()

pdf('pdf/reduced_dim_TSNE_perplexity15.pdf', useDingbats = FALSE)
plotReducedDim(sce_layer, dimred = 'TSNE_perplexity15', colour_by = 'layer_guess') +
    ggplot2::scale_fill_manual(values =  unname(Polychrome::palette36.colors(7)),
        name = 'Layer')
plotReducedDim(sce_layer, dimred = 'TSNE_perplexity15', colour_by = 'sample_name')
plotReducedDim(sce_layer, dimred = 'TSNE_perplexity15', colour_by = 'subject_position')
plotReducedDim(sce_layer, dimred = 'TSNE_perplexity15', colour_by = 'subject')
plotReducedDim(sce_layer, dimred = 'TSNE_perplexity15', colour_by = 'position')
dev.off()

pdf('pdf/reduced_dim_TSNE_perplexity20.pdf', useDingbats = FALSE)
plotReducedDim(sce_layer, dimred = 'TSNE_perplexity20', colour_by = 'layer_guess') +
    ggplot2::scale_fill_manual(values =  unname(Polychrome::palette36.colors(7)),
        name = 'Layer')
plotReducedDim(sce_layer, dimred = 'TSNE_perplexity20', colour_by = 'sample_name')
plotReducedDim(sce_layer, dimred = 'TSNE_perplexity20', colour_by = 'subject_position')
plotReducedDim(sce_layer, dimred = 'TSNE_perplexity20', colour_by = 'subject')
plotReducedDim(sce_layer, dimred = 'TSNE_perplexity20', colour_by = 'position')
dev.off()

pdf('pdf/reduced_dim_UMAP_neighbors15.pdf', useDingbats = FALSE)
plotReducedDim(sce_layer, dimred = 'UMAP_neighbors15', colour_by = 'layer_guess') +
    ggplot2::scale_fill_manual(values =  unname(Polychrome::palette36.colors(7)),
        name = 'Layer')
plotReducedDim(sce_layer, dimred = 'UMAP_neighbors15', colour_by = 'sample_name')
plotReducedDim(sce_layer, dimred = 'UMAP_neighbors15', colour_by = 'subject_position')
plotReducedDim(sce_layer, dimred = 'UMAP_neighbors15', colour_by = 'subject')
plotReducedDim(sce_layer, dimred = 'UMAP_neighbors15', colour_by = 'position')
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

## Save for later
save(sce_layer, file = 'rda/sce_layer.Rdata')
## For mapping back to the original sce object
save(layerIndexes, file = 'rda/layerIndexes.Rdata')
## For subsetting again the genes if necessary
save(ix_mito, selected_genes, file = 'rda/selected_genes.Rdata')


## From here('Analysis', 'convert_sce.R')
sort_clusters <- function(clusters, map_subset = NULL) {
    if (is.null(map_subset)) {
        map_subset <- rep(TRUE, length(clusters))
    }
    map <-
        rank(length(clusters[map_subset]) - table(clusters[map_subset]), ties.method = 'first')
    res <- map[clusters]
    factor(res)
}

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
            gene.names = rowData(sce_layer)$gene_name
        )
    })
    names(res_direc) <- directions
    return(res_direc)
})
names(markers_layer) <- c('any', 'all')
# 2020-01-22 10:46:09 finding markers for p-val any and any direction
# 2020-01-22 10:46:18 finding markers for p-val any and up direction
# 2020-01-22 10:46:27 finding markers for p-val all and any direction
# 2020-01-22 10:46:36 finding markers for p-val all and up direction
save(markers_layer, file = 'rda/markers_layer.Rdata')


## I downloaded the current git version (026edd8eb68fa0479769450473a31795ae50c742)
## to obtain the code for this function
## git clone https://git.bioconductor.org/packages/scran
getMarkerEffects <- function(x, prefix = "logFC", strip = TRUE) {
    regex <- paste0("^", prefix, "\\.")
    i <- grep(regex, colnames(x))
    out <- as.matrix(x[, i])
    
    if (strip) {
        colnames(out) <- sub(regex, "", colnames(out))
    }
    out
}

## Read in the pieces for the gene annotation below
genes_km_raw <-
    read_xlsx(here('Analysis', 'KRM_Layer_Markers.xlsx'))
genes_bm_raw <-
    read_xlsx(here('cortical layer marker gene list_1.xlsx'))


## Build gene annotation data.frame for heatmap
gene_ann <- function(x) {
    m_km <- match(tolower(x), tolower(genes_km_raw$Gene))
    m_bm <- match(tolower(x), tolower(genes_bm_raw$Gene))
    res <-
        data.frame(KM_Zeng = factor(!is.na(m_km), levels = c('FALSE', 'TRUE')),
            BM = factor(!is.na(m_bm), levels = c('FALSE', 'TRUE')))
    rownames(res) <- make.names(x, unique = TRUE)
    return(res)
}
summary(gene_ann(rowData(sce_layer)$gene_name))
#  KM_Zeng          BM
# FALSE:22255   FALSE:22275
# TRUE :   76   TRUE :   56
nrow(genes_km_raw)
# [1] 81
nrow(genes_bm_raw)
# [1] 65
## So 5 and 9 genes are not here to begin with

summary(gene_ann(rowData(sce_original)$gene_name))
#  KM_Zeng          BM
# FALSE:33461   FALSE:33481
# TRUE :   77   TRUE :   57
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


ann_colors <-
    list(
        BM = c(
            'FALSE' = RColorBrewer::brewer.pal(4, 'Dark2')[1],
            'TRUE' = RColorBrewer::brewer.pal(4, 'Dark2')[2]
        ),
        KM_Zeng = c(
            'FALSE' = RColorBrewer::brewer.pal(4, 'Dark2')[3],
            'TRUE' = RColorBrewer::brewer.pal(4, 'Dark2')[4]
        )
    )

## Plotting code
plot_markers_logfc <- function(x, pval.type = c('any', 'all')) {
    lapply(seq_along(x), function(chosen) {
        interesting <- x[[chosen]]
        if (pval.type == 'any') {
            best.set <-
                interesting[interesting$Top <= 6, ] ## for pval.type = 'any'
        } else {
            best.set <- head(interesting, 30) ## for pval.type == 'all'
        }
        logFCs <- getMarkerEffects(best.set)
        print(
            pheatmap(
                logFCs,
                breaks = seq(-3, 3, length.out = 101),
                main = names(x)[chosen],
                # color = colorRampPalette(c("white", "blue"))(100),
                annotation_colors = ann_colors,
                annotation_row = gene_ann(rownames(logFCs)),
                annotation_names_row = TRUE
            )
        )
        return(NULL)
    })
}

sce_layer_symbol <- sce_layer
rownames(sce_layer_symbol) <-
    make.names(rowData(sce_layer)$gene_name, unique = TRUE)

plot_markers_expr <- function(x, pval.type = c('any', 'all')) {
    lapply(seq_along(x), function(chosen) {
        interesting <- x[[chosen]]
        if (pval.type == 'any') {
            best.set <-
                interesting[interesting$Top <= 6, ] ## for pval.type = 'any'
        } else {
            best.set <- head(interesting, 30) ## for pval.type == 'all'
        }
        pheat <- plotHeatmap(
            sce_layer_symbol,
            features = rownames(best.set),
            main = names(x)[chosen],
            colour_columns_by = c('layer_guess', 'subject', 'subject_position', 'sample_name'),
            # color = colorRampPalette(c("white", "blue"))(100),
            # annotation_colors = ann_colors,
            annotation_row = gene_ann(rownames(best.set)),
            annotation_names_row = TRUE
        )
        
        ## Fix the colors for the layers
        layout_num <-
            which(pheat$gtable$layout$name == 'col_annotation')
        layer_names <-
            rownames(pheat$gtable$grobs[[layout_num]]$gp$fill)
        new_layer_cols <-
            Polychrome::palette36.colors(7)[as.integer(sce_layer$layer_guess[match(layer_names, colnames(sce_layer))])]
        names(new_layer_cols) <- layer_names
        pheat$gtable$grobs[[layout_num]]$gp$fill[, 'layer_guess'] <-
            new_layer_cols
        
        ## Now fix the legend
        layout_num <-
            which(pheat$gtable$layout$name == 'annotation_legend')
        children_name <-
            pheat$gtable$grobs[[layout_num]]$childrenOrder['layer_guess r']
        layer_names <-
            names(pheat$gtable$grobs[[layout_num]]$children[[children_name]]$gp$fill)
        new_layer_cols <- Polychrome::palette36.colors(7)
        names(new_layer_cols) <- layer_names
        pheat$gtable$grobs[[layout_num]]$children[[children_name]]$gp$fill <-
            new_layer_cols
        
        ## Print the heatmap
        print(pheat)
        return(NULL)
    })
}

plot_markers_loop <- function(pdf_header, FUN) {
    for (pval in names(markers_layer)) {
        for (direc in names(markers_layer[[1]])) {
            pdf(
                paste0(
                    'pdf/',
                    pdf_header,
                    '_pval_',
                    pval,
                    '_direc_',
                    direc,
                    '.pdf'
                ),
                useDingbats = FALSE,
                height = 14
            )
            FUN(markers_layer[[pval]][[direc]], pval.type = pval)
            dev.off()
        }
    }
}

plot_markers_loop('markers_logFC', plot_markers_logfc)
plot_markers_loop('markers_expr', plot_markers_expr)


## Reproducibility information
print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()

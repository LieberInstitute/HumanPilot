### MNT application of snRNA-seq pilot workflow to Visium DLPFC data
### Initiated MNT 29May2020
### Original SCE file in:
###     /dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Human_DLPFC_Visium_processedData_sce_scran.Rdata
### Working in
###     /dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/spotDissection_MNT/
#####################################################################

library(SingleCellExperiment)
library(EnsDb.Hsapiens.v86)
library(scater)
library(scran)
library(batchelor)
library(DropletUtils)
library(jaffelab)
library(dendextend)
library(dynamicTreeCut)
library(sessioninfo)

### Palette taken from `scater`
tableau10medium = c("#729ECE", "#FF9E4A", "#67BF5C", "#ED665D",
                    "#AD8BC9", "#A8786E", "#ED97CA", "#A2A2A2",
                    "#CDCC5D", "#6DCCDA")
tableau20 = c("#1F77B4", "#AEC7E8", "#FF7F0E", "#FFBB78", "#2CA02C",
              "#98DF8A", "#D62728", "#FF9896", "#9467BD", "#C5B0D5",
              "#8C564B", "#C49C94", "#E377C2", "#F7B6D2", "#7F7F7F",
              "#C7C7C7", "#BCBD22", "#DBDB8D", "#17BECF", "#9EDAE5")

# ===

## Load SCE
load("/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Human_DLPFC_Visium_processedData_sce_scran.Rdata",
     verbose=T)
    # sce, top.hvgs
    length(top.hvgs)  # 1942
  
    # 'top.hvgs' comes from Leo's run in 'sce_scran.R' in the parent dir
        # top.hvgs <- getTopHVGs(dec, prop = 0.1)
        # length(top.hvgs)
        # # [1] 1942

## For reference:
    table(sce$subject_position, sce$replicate)
        #                  1    2
        # Br5292_pos0   4226 4384
        # Br5292_pos300 4789 4634
        # Br5595_pos0   3661 3498
        # Br5595_pos300 4110 4015
        # Br8100_pos0   3639 3673
        # Br8100_pos300 3592 3460

## Make SCE subsets
table(sce$cell_count)
    #    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15
    # 4883 7179 9186 8491 6302 4146 2644 1698 1003  742  453  310  193  149   98   73
    #   16   17   18   19   20   21   22   23   24   25   26   27
    #   44   26   20   17    8    5    2    2    2    2    2    1

sce.singlet <- sce[ ,sce$cell_count==1]
sce.neuropil <- sce[ ,sce$cell_count==0]

# And basically remove all reducedDims
reducedDims(sce.singlet) <- NULL
reducedDims(sce.neuropil) <- NULL


geneVar.singlet <- modelGeneVar(sce.singlet)
chosen.hvgs.singlet <- geneVar.singlet$bio > 0
sum(chosen.hvgs.singlet)  # 15517

geneVar.neuropil <- modelGeneVar(sce.neuropil)
chosen.hvgs.neuropil <- geneVar.neuropil$bio > 0
sum(chosen.hvgs.neuropil) # 14476

    # Check out mean-var trends
    fit.singlet <- metadata(geneVar.singlet)
    plot(fit.singlet$mean, fit.singlet$var, xlab="Mean of log-expression",
         ylab="Variance of log-expression")
    curve(fit.singlet$trend(x), col="dodgerblue", add=TRUE, lwd=2)
    
    fit.neuropil <- metadata(geneVar.neuropil)
    plot(fit.neuropil$mean, fit.neuropil$var, xlab="Mean of log-expression",
         ylab="Variance of log-expression")
    curve(fit.neuropil$trend(x), col="dodgerblue", add=TRUE, lwd=2)
        # These look fine

    
# Go ahead and save these - we'll work with the 'logcounts' already generated
save(sce.singlet, chosen.hvgs.singlet, file="./SCE_singlet-spots_MNT.rda")
save(sce.neuropil, chosen.hvgs.neuropil, file="./SCE_neuropil-spots_MNT.rda")
    
    

### Dimensionality reduction ================================================================
# Run PCA, taking top 100 (instead of default 50 PCs)
set.seed(109)
sce.singlet <- runPCA(sce.singlet, subset_row=chosen.hvgs.singlet, ncomponents=100,
                  BSPARAM=BiocSingular::RandomParam())
set.seed(109)
sce.neuropil <- runPCA(sce.neuropil, subset_row=chosen.hvgs.neuropil, ncomponents=100,
                      BSPARAM=BiocSingular::RandomParam())


## Make t-SNE/UMAP just with default params
# t-SNE
set.seed(109)
sce.singlet <- runTSNE(sce.singlet, dimred="PCA")

set.seed(109)
sce.neuropil <- runTSNE(sce.neuropil, dimred="PCA")


# UMAP
set.seed(109)
sce.singlet <- runUMAP(sce.singlet, dimred="PCA")

set.seed(109)
sce.neuropil <- runUMAP(sce.neuropil, dimred="PCA")


# Save
save(sce.singlet, chosen.hvgs.singlet, file="./SCE_singlet-spots_MNT.rda")
save(sce.neuropil, chosen.hvgs.neuropil, file="./SCE_neuropil-spots_MNT.rda")



### Clustering: Two-step ===========================================
### Step 1: Perform graph-based clustering in this PC space
#         - take k=20 NN to build graph
snn.gr.singlet <- buildSNNGraph(sce.singlet, k=20, use.dimred="PCA")
clusters.k20.singlet <- igraph::cluster_walktrap(snn.gr.singlet)$membership
table(clusters.k20.singlet)
    #    1    2    3    4    5    6    7    8
    #  441  895 2351   43 1408 1945   86   10

# Is sample/position driving any clusters at this level?
table(sce.singlet$subject_position, clusters.k20.singlet)
    #               clusters.k20.singlet
    #                  1    2    3    4    5    6    7    8
    # Br5292_pos0    196  223  580    2  324  126    7    0
    # Br5292_pos300  169  316  533   14  177   33    1    8
    # Br5595_pos0     32  107  112    0   99  525    5    0
    # Br5595_pos300   12  109  134    0  258 1182   41    0
    # Br8100_pos0      1   56  318    8  226    9   12    1
    # Br8100_pos300   31   84  674   19  324   70   20    1

# Assign as 'prelimCluster'
sce.singlet$prelimCluster_MNT <- factor(clusters.k20.singlet)



## Do the same for neuropil
snn.gr.neuropil <- buildSNNGraph(sce.neuropil, k=20, use.dimred="PCA")
clusters.k20.neuropil <- igraph::cluster_walktrap(snn.gr.neuropil)$membership
table(clusters.k20.neuropil)
    #   1    2    3    4    5    6    7    8    9
    #1050   22   85    8  827 1705  688  492    6

# Is sample/position driving any clusters at this level?
table(sce.neuropil$subject_position, clusters.k20.neuropil)
    #                 1   2   3   4   5   6   7   8   9
    # Br5292_pos0   661   0   0   0 218 276 295 154   0
    # Br5292_pos300 281  19  85   7 170 360 232 337   6
    # Br5595_pos0    64   0   0   0  67 129  33   0   0
    # Br5595_pos300  13   0   0   0 245 427  37   1   0
    # Br8100_pos0     4   1   0   1  66 127  39   0   0
    # Br8100_pos300  27   2   0   0  61 386  52   0   0

# Assign as 'prelimCluster'
sce.neuropil$prelimCluster_MNT <- factor(clusters.k20.neuropil)

# Save
save(sce.singlet, chosen.hvgs.singlet, file="./SCE_singlet-spots_MNT.rda")
save(sce.neuropil, chosen.hvgs.neuropil, file="./SCE_neuropil-spots_MNT.rda")

### ** Skipping hierarchical cluster for now because these cluster numbers are small
  #    (approach was conceived bc originally we were working with ~30-200 'prelimClusters',
  #     depending on the setting--region-specific vs pan-brain)

# How do these look in TSNE space?
plotTSNE(sce.singlet, colour_by="prelimCluster_MNT")
plotTSNE(sce.neuropil, colour_by="prelimCluster_MNT")
    # Neither of these look great


head(attr(reducedDim(sce.singlet, "PCA"), "percentVar"), n=10)
    # [1] 1.15381985 0.56135136 0.41521433 0.38312903 0.28329346 0.22935132
    # [7] 0.17669649 0.15619886 0.14740650 0.14509740
sum(attr(reducedDim(sce.singlet, "PCA"), "percentVar")) # so only 11.1% var across 100 PCs


head(attr(reducedDim(sce.neuropil, "PCA"), "percentVar"), n=10)
    # [1] 1.2675013 0.6557071 0.5020080 0.4714249 0.3193334 0.2207257 0.2088314
    # [8] 0.1806465 0.1609374 0.1496319 
sum(attr(reducedDim(sce.neuropil, "PCA"), "percentVar"))  # only 12.9% var across 100 PCs

    ## What were these look like in the original SCE?
    head(attr(reducedDim(sce, "PCA"), "percentVar"), n=10)
        # [1] 3.7490227 1.8043903 1.2798964 0.8183388 0.6291475 0.4590034 0.4115223
        # [8] 0.3765352 0.3022602 0.2916951
    
    sum(attr(reducedDim(sce, "PCA"), "percentVar")) # 17.5%   - so not too much better

        # MNT remarks: in snRNA-seq datasets, usually PC1 alone is ~15-20% var and the top
        #              100 PCs usually account for ~40-50% total var


    
## Let's go ahead and plot some expression of these clusters:
    
# To make more easily accessible...
rowData(sce.singlet)$gene_name_unique <- uniquifyFeatureNames(rowData(sce.singlet)$gene_id, rowData(sce.singlet)$gene_name)
rownames(sce.singlet) <- rowData(sce.singlet)$gene_name_unique
rowData(sce.neuropil)$gene_name_unique <- uniquifyFeatureNames(rowData(sce.neuropil)$gene_id, rowData(sce.neuropil)$gene_name)
rownames(sce.neuropil) <- rowData(sce.neuropil)$gene_name_unique
# Marker genes as used in snRNA-seq workflow
markers.mathys.custom = list(
  'neurons' = c('SYT1', 'SNAP25', 'GRIN1'),
  'excitatory_neuron' = c('CAMK2A', 'NRGN','SLC17A7', 'SLC17A6', 'SLC17A8'),
  'inhibitory_neuron' = c('GAD1', 'GAD2', 'SLC32A1'),
  'oligodendrocyte' = c('MBP', 'MOBP', 'PLP1'),
  'oligodendrocyte_precursor' = c('PDGFRA', 'VCAN', 'CSPG4'),
  'microglia' = c('CD74', 'CSF1R', 'C3'),
  'astrocytes' = c('GFAP', 'TNC', 'AQP4', 'SLC1A2'),
  'endothelial' = c('CLDN5', 'FLT1', 'VTN')
)

# Singlet spots
pdf("./pdfs/broadMarker-logExprs_graphClusters_singlet-spots_MNT.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.singlet, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
                   x="prelimCluster_MNT", colour_by="prelimCluster_MNT", point_alpha=0.5, point_size=.7,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium[1:8], length(markers.mathys.custom[[i]]))) +  
      ggtitle(label=paste0(names(markers.mathys.custom)[i], " broad markers"))
  )
}
dev.off()

    # For reviewers 03Sep2020:
    library(gridExtra)
    snap25 <- plotExpression(sce.singlet, exprs_values = "logcounts", features="SNAP25",
                             x="prelimCluster_MNT", colour_by="prelimCluster_MNT", point_alpha=0.5, point_size=.7,
                             add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                          geom = "crossbar", width = 0.3,
                                                          colour=tableau10medium[1:8]) +
              xlab("")
    
    mbp <- plotExpression(sce.singlet, exprs_values = "logcounts", features="MBP",
                             x="prelimCluster_MNT", colour_by="prelimCluster_MNT", point_alpha=0.5, point_size=.7,
                             add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                          geom = "crossbar", width = 0.3,
                                                          colour=tableau10medium[1:8]) +
           xlab("")
    
    pdf("./pdfs/SNAP25-MBP_singlet-spots-clustered_MNT.pdf", height=6, width=4)
    grid.arrange(snap25, mbp, ncol=1)
    dev.off()



# Neuropil spots
pdf("./pdfs/broadMarker-logExprs_graphClusters_neuropil-spots_MNT.pdf", height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpression(sce.neuropil, exprs_values = "logcounts", features=c(markers.mathys.custom[[i]]),
                   x="prelimCluster_MNT", colour_by="prelimCluster_MNT", point_alpha=0.5, point_size=.7,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium[1:9], length(markers.mathys.custom[[i]]))) +  
      ggtitle(label=paste0(names(markers.mathys.custom)[i], " broad markers"))
  )
}
dev.off()

## Also some requested markers by ST group ===
    g_pre = toupper(c("Bsn", "Rims1", "Rims2", "Rims3", "Stx6", "Snap25","Rapgef4"))
    g_post = toupper(c("Arc", "Bdnf", "Grin1", "SLC17A7", "Gria1"))
    g_other = toupper(c("Dlg4", "Grin2a", "Limk1", "Fmr1"))
    g_ecm = toupper(c("HAPLN1", "ACAN", "BCAN", "TNR"))
    g_my = toupper(c("Mobp", "Mbp", "Aqp4"))
    ix_mito <- grep("^MT-", rowData(sce.singlet)$gene_name)
    g_mito = rowData(sce.singlet)$gene_name[ix_mito]
    
    g = list("g_pre"=g_pre,
             "g_post"=g_post,
             "g_other"=g_other,
             "g_ecm"=g_ecm,
             "g_my"=g_my,
             "g_mito"=g_mito)
    table(unlist(g) %in% rownames(sce.singlet)) # all there
    
    # Singlet
    pdf("./pdfs/selectMarkers-AnJa-logExprs_graphClusters_singlet-spots_MNT.pdf", height=5, width=8)
    for(i in 1:length(g)){
      print(
        plotExpression(sce.singlet, exprs_values = "logcounts", features=c(g[[i]]),
                       x="prelimCluster_MNT", colour_by="prelimCluster_MNT", point_alpha=0.5, point_size=.7,
                       ncol=3, add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                    geom = "crossbar", width = 0.3,
                                                    colour=rep(tableau10medium[1:8], length(g[[i]]))) +  
          ggtitle(label=paste0(names(g)[i], ": select markers"))
      )
    }
    dev.off()
    
    # Neuropil
    pdf("./pdfs/selectMarkers-AnJa-logExprs_graphClusters_neuropil-spots_MNT.pdf", height=5, width=8)
    for(i in 1:length(g)){
      print(
        plotExpression(sce.neuropil, exprs_values = "logcounts", features=c(g[[i]]),
                       x="prelimCluster_MNT", colour_by="prelimCluster_MNT", point_alpha=0.5, point_size=.7,
                       ncol=3, add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                    geom = "crossbar", width = 0.3,
                                                    colour=rep(tableau10medium[1:9], length(g[[i]]))) +  
          ggtitle(label=paste0(names(g)[i], " broad markers"))
      )
    }
    dev.off()


## Let's plot some reducedDims too ===
# First set the 'sample_name' to character
sce.singlet$sample_name <- as.character(sce.singlet$sample_name)
sce.neuropil$sample_name <- as.character(sce.neuropil$sample_name)


pdf("./pdfs/reducedDims-with-MNT-prelimClusters_singlet-spots.pdf")
plotReducedDim(sce.singlet, dimred="PCA", ncomponents=5, colour_by="prelimCluster_MNT", point_alpha=0.5)
plotTSNE(sce.singlet, colour_by="subject_position", point_size=3.5, point_alpha=0.5)
plotTSNE(sce.singlet, colour_by="sample_name", point_size=3.5, point_alpha=0.5)
plotTSNE(sce.singlet, colour_by="prelimCluster_MNT", point_size=3.5, point_alpha=0.5)
dev.off()

pdf("./pdfs/reducedDims-with-MNT-prelimClusters_neuropil-spots.pdf")
plotReducedDim(sce.neuropil, dimred="PCA", ncomponents=5, colour_by="prelimCluster_MNT", point_alpha=0.5)
plotTSNE(sce.neuropil, colour_by="subject_position", point_size=3.5, point_alpha=0.5)
plotTSNE(sce.neuropil, colour_by="sample_name", point_size=3.5, point_alpha=0.5)
plotTSNE(sce.neuropil, colour_by="prelimCluster_MNT", point_size=3.5, point_alpha=0.5)
dev.off()


# Save
save(sce.singlet, chosen.hvgs.singlet, file="./SCE_singlet-spots_MNT.rda")
save(sce.neuropil, chosen.hvgs.neuropil, file="./SCE_neuropil-spots_MNT.rda")





## ** Exploration: Should re-generate log-normalized counts prior to going into these analyses?
cell_count_idx <- splitit(sce$cell_count)
sapply(cell_count_idx, function(x){quantile(sce[ ,x]$sum_umi)})
    #           0     1     2       3       4        5        6        7       8
    # 0%      17.0    58    58   126.0    74.0   151.00   140.00   298.00   264.0
    # 25%   1338.0  1891  2137  2220.5  2301.5  2333.25  2248.25  2253.50  2230.5
    # 50%   2058.0  2762  3097  3260.0  3383.0  3449.00  3301.00  3273.50  3282.0
    # 75%   3006.5  3923  4474  4763.0  4824.5  4962.50  4765.25  4743.75  4669.5
    # 100% 11163.0 14904 18266 17295.0 17675.0 20059.00 20600.00 17438.00 15273.0
    #            9    10       11   12   13      14   15     16      17     18   19
    # 0%     209.0   215   374.00  483  675  762.00  663  768.0  641.00  709.0 1129
    # 25%   2265.0  2216  2227.25 1923 2100 2077.75 1550 1892.0 2279.25 1722.5 1736
    # 50%   3258.5  3377  3194.50 2832 2981 3115.50 2707 3128.0 3093.50 2175.0 2238
    # 75%   4461.5  4507  4598.50 3987 3905 4327.00 3515 3641.5 3960.25 3268.0 2779
    # 100% 16151.0 11636 12834.00 8712 8614 8119.00 7682 6638.0 5680.00 4047.0 4198
    #          20   21   22      23     24      25      26  27
    # 0%    781.0 1792 3205 1029.00 1277.0 3091.00 2351.00 840
    # 25%  2153.5 2230 3459 1566.75 1521.5 3192.75 2495.25 840
    # 50%  2766.5 2956 3713 2104.50 1766.0 3294.50 2639.50 840    - observation: the spot library sizes
    # 75%  2995.5 3274 3967 2642.25 2010.5 3396.25 2783.75 840      generally increase with more cells/spot
    # 100% 3925.0 4066 4221 3180.00 2255.0 3498.00 2928.00 840      (this is expected)

sapply(cell_count_idx, function(x){round(quantile(sizeFactors(sce)[x]),3)})
    #         0     1     2     3     4     5     6     7     8     9    10    11
    # 0%   0.001 0.007 0.005 0.020 0.019 0.043 0.032 0.070 0.058 0.064 0.059 0.104
    # 25%  0.345 0.506 0.581 0.612 0.644 0.653 0.630 0.631 0.624 0.635 0.598 0.620
    # 50%  0.550 0.759 0.873 0.922 0.972 0.999 0.952 0.945 0.946 0.927 0.943 0.902
    # 75%  0.840 1.106 1.283 1.378 1.410 1.482 1.424 1.398 1.365 1.314 1.309 1.298
    # 100% 3.395 4.673 5.949 5.764 5.955 6.760 6.913 6.023 4.926 5.060 3.691 4.228
    # 12    13    14    15    16    17    18    19    20    21    22    23
    # 0%   0.114 0.170 0.242 0.175 0.215 0.177 0.210 0.327 0.224 0.498 0.936 0.298
    # 25%  0.551 0.599 0.604 0.457 0.572 0.640 0.443 0.500 0.628 0.731 1.027 0.452
    # 50%  0.807 0.844 0.870 0.775 0.851 0.910 0.610 0.679 0.788 0.903 1.119 0.605
    # 75%  1.124 1.103 1.196 1.046 1.024 1.138 0.895 0.819 0.890 0.921 1.210 0.759
    # 100% 2.847 2.526 2.631 2.323 2.076 1.710 1.257 1.274 1.000 1.156 1.302 0.912
    # 24    25    26    27
    # 0%   0.363 0.940 0.733 0.233
    # 25%  0.434 0.948 0.779 0.233
    # 50%  0.504 0.957 0.824 0.233      - oh so for where cell_count == 0 & 1, the logcounts
    # 75%  0.575 0.966 0.869 0.233        are 'inflated', so to speak, so this is fine. (Would have
    # 100% 0.646 0.975 0.914 0.233        wanted to re-compute if the normalizing sizeFactors were > 1)

# How do de-novo LSFs correlate with existing?
cor(sizeFactors(sce)[cell_count_idx[["0"]]], librarySizeFactors(sce.neuropil))  # [1] 0.9917861
cor(sizeFactors(sce)[cell_count_idx[["1"]]], librarySizeFactors(sce.singlet)) # [1] 0.9933599

pdf("pdfs/LSFcomparisoin_deNovo-vs-exististing_MNT.pdf")
plot(sizeFactors(sce)[cell_count_idx[["0"]]], librarySizeFactors(sce.neuropil),
     xlab="Size factors computed across all Visum spots (cell_count=0)", ylab="de-novo LSFs",
     main="Library Size Factor comparison: neuropil \n (scaled to mean SFs across spots)")
plot(sizeFactors(sce)[cell_count_idx[["1"]]], librarySizeFactors(sce.singlet),
     xlab="Size factors computed across all Visum spots (cell_count=1)", ylab="de-novo LSFs",
     main="Library Size Factor comparison: singlet spots \n (scaled to mean SFs across spots)")  
dev.off()



### Marker tests ================================================================================
    ## Run ANOVA real quick 
    library(edgeR)
    library(doMC)
    registerDoMC(cores=4)
    
    # Take non-0-count genes only
    mat = assays(sce.singlet)$logcounts[rowSums(assay(sce.singlet, "counts"))!=0, ] # 22355 genes
    
    # Do regression
    varCompAnalysis.subject = foreach(i = 1:nrow(mat)) %dopar% {
      if(i %% 1000 == 0) cat("...")
      fit = lm(as.numeric(mat[i,]) ~ subject + position,
               data=colData(sce.singlet))
      full = anova(fit)
      fullSS = full$"Sum Sq"
      signif(cbind(full, PctExp = fullSS/sum(fullSS)*100), 3)
    }
    names(varCompAnalysis.subject) = rownames(mat)
    
    # or each individual sample
    varCompAnalysis.12sxns = foreach(i = 1:nrow(mat)) %dopar% {
      if(i %% 1000 == 0) cat("...")
      fit = lm(as.numeric(mat[i,]) ~ sample_name,
               data=colData(sce.singlet))
      full = anova(fit)
      fullSS = full$"Sum Sq"
      signif(cbind(full, PctExp = fullSS/sum(fullSS)*100), 3)
    }
    names(varCompAnalysis.12sxns) = rownames(mat)
    
    ## make boxplot
    varExpl.subject = t(sapply(varCompAnalysis.subject, function(x) x[,"PctExp"]))
    colnames(varExpl.subject) = rownames(varCompAnalysis.subject[[1]])
    
        round(quantile(varExpl.subject[ ,"subject"], probs=seq(0.1,1,by=0.1)),3)
        #    10%    20%    30%    40%    50%    60%    70%    80%    90%   100%
        #  0.007  0.015  0.023  0.029  0.040  0.054  0.079  0.121  0.225 73.200
    
    varExpl.12sxns = t(sapply(varCompAnalysis.12sxns, function(x) x[,"PctExp"]))
    colnames(varExpl.12sxns) = rownames(varCompAnalysis.12sxns[[1]])
        round(quantile(varExpl.12sxns[ ,"sample_name"], probs=seq(0.1,1,by=0.1)),3)
        #   10%    20%    30%    40%    50%    60%    70%    80%    90%   100%
        # 0.098  0.122  0.141  0.160  0.186  0.218  0.256  0.319  0.450 73.800
        
        
    pdf("pdfs/anova_singletSpots-bySubject-position_or_section-n12_MNT.pdf")
    boxplot(varExpl.subject, main="ANOVA on single-cell spots: `lm( ~ subject + position)`",
            ylab="Percent Var explained (%))")
    boxplot(varExpl.12sxns, main="ANOVA on single-cell spots: lm( ~ sample_name) (n=12)",
            ylab="Percent Var explained (%))")
    dev.off()
    
    save(varCompAnalysis.subject, varCompAnalysis.12sxns, file="./anova_output_subject-pos_or_section-n12_MNT.rda")


    
    
## Let's go ahead and get markers, modeling across sample names/sections for singlet spots ============
    
## Pairwise t-tests ===
mod <- with(colData(sce.singlet), model.matrix(~ sample_name))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`
    ## If modeled with 0 intercept (and not dropping - i.e. mod is of 12 columns per 12 sxns),
     # get: "Error in .ranksafe_qr(full.design) : design matrix is not of full rank

# Drop all 0 genes
sce.singlet <- sce.singlet[!rowSums(assay(sce.singlet,"counts"))==0, ]  # 22355 x 7179

# Run pairwise t-tests
markers.singlet.t.pw <- findMarkers(sce.singlet, groups=sce.singlet$prelimCluster_MNT,
                                    assay.type="logcounts", design=mod, test="t",
                                    direction="up", pval.type="all", full.stats=T)

sapply(markers.singlet.t.pw, function(x){table(x$FDR<0.05)})
    # Not many - a few to a few dozen - otherwise 394 for cluster 8 (which is entirely 'WM')


## Cluster-vs-all-others iteration ===
markers.singlet.t.1vAll <- list()
for(i in levels(sce.singlet$prelimCluster_MNT)){
  # Make temporary contrast
  sce.singlet$contrast <- ifelse(sce.singlet$prelimCluster_MNT==i, 1, 0)
  # Test cluster vs. all
  markers.singlet.t.1vAll[[i]] <- findMarkers(sce.singlet, groups=sce.singlet$contrast,
                                          assay.type="logcounts", design=mod, test="t",
                                          direction="up", pval.type="all", full.stats=T)
}

## Then, temp set of stats to get the standardized logFC
temp.1vAll <- list()
for(i in levels(sce.singlet$prelimCluster_MNT)){
  # Make temporary contrast
  sce.singlet$contrast <- ifelse(sce.singlet$prelimCluster_MNT==i, 1, 0)
  # Test cluster vs. all
  temp.1vAll[[i]] <- findMarkers(sce.singlet, groups=sce.singlet$contrast,
                                 assay.type="logcounts", design=mod, test="t",
                                 std.lfc=TRUE,
                                 direction="up", pval.type="all", full.stats=T)
}


    ## For some reason all the results are in the second List entry (first is always empty)

# Replace that empty slot with the entry with the actul stats
markers.singlet.t.1vAll <- lapply(markers.singlet.t.1vAll, function(x){ x[[2]] })
# Same for that with std.lfc
temp.1vAll <- lapply(temp.1vAll, function(x){ x[[2]] })

# Now just pull from the 'stats.0' DataFrame column
markers.singlet.t.1vAll <- lapply(markers.singlet.t.1vAll, function(x){ x$stats.0 })
temp.1vAll <- lapply(temp.1vAll, function(x){ x$stats.0 })

# Re-name std.lfc column and add to the first result
for(i in names(temp.1vAll)){
  colnames(temp.1vAll[[i]])[1] <- "std.logFC"
  markers.singlet.t.1vAll[[i]] <- cbind(markers.singlet.t.1vAll[[i]], temp.1vAll[[i]]$std.logFC)
  # Oh the colname is kept weird
  colnames(markers.singlet.t.1vAll[[i]])[4] <- "std.logFC"
  # Then re-organize
  markers.singlet.t.1vAll[[i]] <- markers.singlet.t.1vAll[[i]][ ,c("logFC","std.logFC","log.p.value","log.FDR")]
}
    

## Let's save this along with the previous pairwise results
save(markers.singlet.t.1vAll, markers.singlet.t.pw,
     file="markers-stats_singlet-spot-MNT-clusters_findMarkers.rda")




## Print these to pngs
markerList.t.1vAll <- lapply(markers.singlet.t.1vAll, function(x){
  rownames(x)[x[ ,"log.FDR"] < log10(0.000001)]
  }
)
lengths(markerList.t.1vAll)
    #   1    2    3    4    5    6    7    8
    # 356  191 1243  595  179  787  317  494

genes.top40.t.1vAll <- lapply(markerList.t.1vAll, function(x){head(x, n=40)})

for(i in names(genes.top40.t.1vAll)){
  png(paste0("pdfs/markers_singletSpots_1vALL_top40markers_prelimCluster.",i,"_logExprs.png"), height=1900, width=1200)
  print(
    plotExpression(sce.singlet, exprs_values = "logcounts", features=genes.top40.t.1vAll[[i]],
                   x="prelimCluster_MNT", colour_by="prelimCluster_MNT", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium[1:8], length(genes.top40.t.1vAll[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0("Singlet-cell spot prelimCluster_MNT:",i, " top 40 markers: t-tests, cluster-vs-all-others"))
  )
  dev.off()
}





## And for those with any significant pairwise results, print those [up to 40] markers too ===
# Print these to pngs
markerList.t.pw <- lapply(markers.singlet.t.pw, function(x){
  rownames(x)[x$FDR < 0.05]
  }
)

## First write these top 40 lists to a csv
genes.top40.t.pw <- lapply(markerList.t.pw, function(x){head(x, n=40)})

# PW result for those doesn't have 40 markers:
for(i in 1:length(genes.top40.t.pw)){
  genes.top40.t.pw[[i]] <- if(length(genes.top40.t.pw[[i]]) < 40){
    c(genes.top40.t.pw[[i]], rep("",40-length(genes.top40.t.pw[[i]])))
  } else {
    genes.top40.t.pw[[i]]
  }
}

names(genes.top40.t.pw) <- paste0(names(genes.top40.t.pw),"_pw")
names(markerList.t.1vAll) <- paste0(names(markerList.t.1vAll),"_1vAll")

top40genes <- cbind(sapply(genes.top40.t.pw, cbind),
                    sapply(markerList.t.1vAll, function(y) head(y, n=40)))
top40genes <- top40genes[ ,sort(colnames(top40genes))]

write.csv(top40genes, file="./top40genesLists_singletSpot-prelimCluster-MNT_Jun2020.csv",
          row.names=FALSE)



## Will have to print these manually bc such drastically different numbers of significantly-enriched genes:
genes.top40.t.pw <- lapply(markerList.t.pw, function(x){head(x, n=40)})
lengths(genes.top40.t.pw)
    #  1  2  3  4  5  6  7  8
    # 15  0  1 37  0  7  0 40

# 4 & 8:
for(i in c(4,8)){
  png(paste0("pdfs/markers_singletSpots_pairwise_top40markers_prelimCluster.",i,"_logExprs.png"), height=1900, width=1200)
  print(
    plotExpression(sce.singlet, exprs_values = "logcounts", features=genes.top40.t.pw[[i]],
                   x="prelimCluster_MNT", colour_by="prelimCluster_MNT", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium[1:8], length(genes.top40.t.pw[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0("Singlet-cell spot prelimCluster_MNT:",i, " top 40 markers: t-tests, pairwise"))
  )
  dev.off()
}
# 1, 3, & 6
for(i in c(1,3,6)){
  png(paste0("pdfs/markers_singletSpots_pairwise_top40markers_prelimCluster.",i,"_logExprs.png"), height=475, width=1200)
  print(
    plotExpression(sce.singlet, exprs_values = "logcounts", features=genes.top40.t.pw[[i]],
                   x="prelimCluster_MNT", colour_by="prelimCluster_MNT", point_alpha=0.5, point_size=.7, ncol=5,
                   add_legend=F) + stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                geom = "crossbar", width = 0.3,
                                                colour=rep(tableau10medium[1:8], length(genes.top40.t.pw[[i]]))) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(size = 25)) +  
      ggtitle(label=paste0("Singlet-cell spot prelimCluster_MNT:",i, " top 40 markers: t-tests, pairwise"))
  )
  dev.off()
}




### MNT: Correlation to DLPFC snRNA-seq subclusters =========
load("./SCE_singlet-spots_MNT.rda")
    # markers.singlet.t.pw, markers.singlet.t.1vAll
    rm(markers.singlet.t.pw)

load("./SCE_singlet-spots_MNT.rda")
    # sce.singlet, chosen.hvgs.singlet

## Calculate and add t-statistic (= std.logFC * sqrt(N)) for mouse clusters
#      and fix row order to the first entry "1"
fixTo <- rownames(markers.singlet.t.1vAll[[1]])
for(x in names(markers.singlet.t.1vAll)){
  markers.singlet.t.1vAll[[x]]$t.stat <- markers.singlet.t.1vAll[[x]]$std.logFC * sqrt(ncol(sce.singlet))
  markers.singlet.t.1vAll[[x]] <- markers.singlet.t.1vAll[[x]][fixTo, ]
}

# Pull out the t's
ts.singlet <- sapply(markers.singlet.t.1vAll, function(x){x$t.stat})
rownames(ts.singlet) <- fixTo


## Bring in human stats; create t's ===
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/markers-stats_DLPFC_n2_findMarkers-SN-LEVEL_MNTApr2020.rda", verbose=T)
    # markers.dlpfc.t.1vAll, markers.dlpfc.t.design
    rm(markers.dlpfc.t.design)

# Need to add t's with N nuclei used in constrasts
load("/dcl01/lieber/ajaffe/Matt/MNT_thesis/snRNAseq/10x_pilot_FINAL/rdas/regionSpecific_DLPFC-n2_SCE_cellTypesSplit-fromST_Apr2020.rda", verbose=T)
    # sce.dlpfc.st, clusterRefTab.dlpfc, chosen.hvgs.dlpfc, ref.sampleInfo
    rm(clusterRefTab.dlpfc, chosen.hvgs.dlpfc, ref.sampleInfo)

# First drop "Ambig.lowNtrxts" (168 nuclei)
sce.dlpfc.st <- sce.dlpfc.st[ ,sce.dlpfc.st$cellType.split != "Ambig.lowNtrxts"]
sce.dlpfc.st$cellType.split <- droplevels(sce.dlpfc.st$cellType.split)

## As above, calculate and add t-statistic (= std.logFC * sqrt(N)) from contrasts
#      and fix row order to the first entry "Astro"
fixTo <- rownames(markers.dlpfc.t.1vAll[["Astro"]])
for(s in names(markers.dlpfc.t.1vAll)){
  markers.dlpfc.t.1vAll[[s]]$t.stat <- markers.dlpfc.t.1vAll[[s]]$std.logFC * sqrt(ncol(sce.dlpfc.st))
  markers.dlpfc.t.1vAll[[s]] <- markers.dlpfc.t.1vAll[[s]][fixTo, ]
}
# Pull out the t's
ts.dlpfc.sn <- sapply(markers.dlpfc.t.1vAll, function(x){x$t.stat})
rownames(ts.dlpfc.sn) <- fixTo


# Take intersecting and match rows
table(rownames(ts.dlpfc.sn) %in% rownames(ts.singlet))  # 21,117
intersect.dlpfc <- intersect(rownames(ts.dlpfc.sn), rownames(ts.singlet))

ts.dlpfc.sn <- ts.dlpfc.sn[intersect.dlpfc, ]
ts.singlet <- ts.singlet[intersect.dlpfc, ]

cor_ts_snVis <- cor(ts.dlpfc.sn, ts.singlet)

signif(cor_ts_snVis, 2)
    #                     1      2      3      4       5      6      7     8
    # Oligo           0.120  0.160 -0.240  0.400  0.0210 -0.130  0.310  0.37
    # Astro           0.270  0.180 -0.170 -0.110 -0.0650 -0.093 -0.050 -0.13
    # Inhib.4        -0.240 -0.250  0.280 -0.290  0.0240  0.210 -0.250 -0.24
    # Excit.L4:5     -0.360 -0.380  0.390 -0.330  0.0200  0.320 -0.320 -0.27
    # Micro           0.100  0.140 -0.110 -0.041  0.0140 -0.090  0.024 -0.06
    # Inhib.6        -0.240 -0.240  0.290 -0.290 -0.0032  0.190 -0.260 -0.25
    # OPC             0.039  0.014  0.024 -0.120 -0.0260 -0.014 -0.070 -0.12
    # Excit.L2:3     -0.310 -0.320  0.420 -0.310 -0.0220  0.180 -0.280 -0.25
    # Excit.ambig    -0.340 -0.350  0.450 -0.330 -0.0130  0.200 -0.300 -0.27
    # Excit.L3:4     -0.300 -0.320  0.370 -0.310 -0.0050  0.240 -0.280 -0.25
    # Excit.L5:6     -0.280 -0.290  0.330 -0.290  0.0230  0.210 -0.230 -0.25
    # Excit.L6.broad -0.330 -0.340  0.390 -0.320  0.0170  0.240 -0.280 -0.27
    # Inhib.5        -0.230 -0.230  0.280 -0.280  0.0094  0.180 -0.240 -0.24
    # Excit.L5       -0.220 -0.240  0.240 -0.240  0.0280  0.200 -0.210 -0.21
    # Inhib.1        -0.170 -0.180  0.210 -0.200 -0.0071  0.160 -0.190 -0.17
    # Inhib.2        -0.200 -0.210  0.220 -0.210  0.0220  0.160 -0.190 -0.18
    # Inhib.3        -0.190 -0.210  0.230 -0.230  0.0130  0.180 -0.220 -0.19



## Print this to heatmap:
theSeq.all = seq(-.5, .5, by = 0.01)
my.col.all <- colorRampPalette(brewer.pal(7, "BrBG"))(length(theSeq.all)-1)

pdf("pdfs/cor-t-heatmap_Vis-singlet-cell-clusters_to_DLPFC-snRNA-seq.pdf")
pheatmap(cor_ts_snVis,
         color=my.col.all,
         cluster_rows=F, #cluster_cols=F,
         angle_col=90,
         breaks=theSeq.all,
         fontsize=10, fontsize_row=15, fontsize_col=18,
         display_numbers=T, number_format="%.2f", fontsize_number=8,
         main="                       Correlation of singlet-cell cluster t's to DLPFC snRNA-seq subclusters")
dev.off()




print('Reproducibility information:')
Sys.time()
proc.time()
options(width = 120)
session_info()
## =====================================
    # [1] "Reproducibility information:"
    # > Sys.time()
    # [1] "2020-05-30 00:58:46 EDT"
    # > proc.time()
    # user   system  elapsed
    # 1407.251   24.439 8130.455
    # > options(width = 120)
    # > session_info()
    # * 0.32.0   2019-10-29 [2] Bioconductor
    # BiocNeighbors          1.4.2    2020-02-29 [1] Bioconductor
    # BiocParallel         * 1.20.1   2019-12-21 [2] Bioconductor
    # BiocSingular           1.2.1    2019-12-23 [1] Bioconductor
    # biomaRt                2.42.0   2019-10-29 [2] Bioconductor
    # Biostrings             2.54.0   2019-10-29 [2] Bioconductor
    # bit                    1.1-15.2 2020-02-10 [2] CRAN (R 3.6.1)
    # bit64                  0.9-7    2017-05-08 [2] CRAN (R 3.6.1)
    # bitops                 1.0-6    2013-08-17 [2] CRAN (R 3.6.1)
    # blob                   1.2.1    2020-01-20 [2] CRAN (R 3.6.1)
    # cli                    2.0.1    2020-01-08 [2] CRAN (R 3.6.1)
    # codetools              0.2-16   2018-12-24 [3] CRAN (R 3.6.1)
    # colorspace             1.4-1    2019-03-18 [2] CRAN (R 3.6.1)
    # cowplot                1.0.0    2019-07-11 [1] CRAN (R 3.6.1)
    # crayon                 1.3.4    2017-09-16 [2] CRAN (R 3.6.1)
    # curl                   4.3      2019-12-02 [2] CRAN (R 3.6.1)
    # DBI                    1.1.0    2019-12-15 [2] CRAN (R 3.6.1)
    # dbplyr                 1.4.2    2019-06-17 [2] CRAN (R 3.6.1)
    # DelayedArray         * 0.12.2   2020-01-06 [2] Bioconductor
    # DelayedMatrixStats     1.8.0    2019-10-29 [2] Bioconductor
    # dendextend           * 1.13.3   2020-02-08 [2] CRAN (R 3.6.1)
    # digest                 0.6.24   2020-02-12 [1] CRAN (R 3.6.1)
    # dplyr                  0.8.4    2020-01-31 [2] CRAN (R 3.6.1)
    # dqrng                  0.2.1    2019-05-17 [1] CRAN (R 3.6.1)
    # DropletUtils         * 1.6.1    2019-10-30 [1] Bioconductor
    # dynamicTreeCut       * 1.63-1   2016-03-11 [1] CRAN (R 3.6.1)
    # edgeR                  3.28.0   2019-10-29 [2] Bioconductor
    # EnsDb.Hsapiens.v86   * 2.99.0   2019-10-23 [1] Bioconductor
    # ensembldb            * 2.10.2   2019-11-20 [2] Bioconductor
    # fansi                  0.4.1    2020-01-08 [2] CRAN (R 3.6.1)
    # farver                 2.0.3    2020-01-16 [2] CRAN (R 3.6.1)
    # GenomeInfoDb         * 1.22.0   2019-10-29 [2] Bioconductor
    # GenomeInfoDbData       1.2.2    2019-10-28 [2] Bioconductor
    # GenomicAlignments      1.22.1   2019-11-12 [2] Bioconductor
    # GenomicFeatures      * 1.38.2   2020-02-15 [2] Bioconductor
    # GenomicRanges        * 1.38.0   2019-10-29 [2] Bioconductor
    # ggbeeswarm             0.6.0    2017-08-07 [1] CRAN (R 3.6.1)
    # ggplot2              * 3.2.1    2019-08-10 [2] CRAN (R 3.6.1)
    # glue                   1.3.1    2019-03-12 [2] CRAN (R 3.6.1)
    # googledrive            1.0.1    2020-05-05 [1] CRAN (R 3.6.1)
    # gridExtra              2.3      2017-09-09 [2] CRAN (R 3.6.1)
    # gtable                 0.3.0    2019-03-25 [2] CRAN (R 3.6.1)
    # HDF5Array              1.14.2   2020-01-31 [2] Bioconductor
    # hms                    0.5.3    2020-01-08 [2] CRAN (R 3.6.1)
    # httr                   1.4.1    2019-08-05 [2] CRAN (R 3.6.1)
    # igraph                 1.2.4.2  2019-11-27 [2] CRAN (R 3.6.1)
    # IRanges              * 2.20.2   2020-01-13 [2] Bioconductor
    # irlba                  2.3.3    2019-02-05 [1] CRAN (R 3.6.1)
    # jaffelab             * 0.99.29  2019-10-23 [1] Github (LieberInstitute/jaffelab@a7d87cb)
    # labeling               0.3      2014-08-23 [2] CRAN (R 3.6.1)
    # lattice                0.20-38  2018-11-04 [3] CRAN (R 3.6.1)
    # lazyeval               0.2.2    2019-03-15 [2] CRAN (R 3.6.1)
    # lifecycle              0.1.0    2019-08-01 [2] CRAN (R 3.6.1)
    # limma                  3.42.2   2020-02-03 [2] Bioconductor
    # locfit                 1.5-9.1  2013-04-20 [2] CRAN (R 3.6.1)
    # magrittr               1.5      2014-11-22 [2] CRAN (R 3.6.1)
    # Matrix                 1.2-17   2019-03-22 [3] CRAN (R 3.6.1)
    # matrixStats          * 0.55.0   2019-09-07 [2] CRAN (R 3.6.1)
    # memoise                1.1.0    2017-04-21 [2] CRAN (R 3.6.1)
    # munsell                0.5.0    2018-06-12 [2] CRAN (R 3.6.1)
    # openssl                1.4.1    2019-07-18 [2] CRAN (R 3.6.1)
    # pillar                 1.4.3    2019-12-20 [2] CRAN (R 3.6.1)
    # pkgconfig              2.0.3    2019-09-22 [2] CRAN (R 3.6.1)
    # plyr                   1.8.5    2019-12-10 [2] CRAN (R 3.6.1)
    # prettyunits            1.1.1    2020-01-24 [2] CRAN (R 3.6.1)
    # progress               1.2.2    2019-05-16 [2] CRAN (R 3.6.1)
    # ProtGenerics           1.18.0   2019-10-29 [2] Bioconductor
    # purrr                  0.3.4    2020-04-17 [1] CRAN (R 3.6.1)
    # R.methodsS3            1.8.0    2020-02-14 [2] CRAN (R 3.6.1)
    # R.oo                   1.23.0   2019-11-03 [2] CRAN (R 3.6.1)
    # R.utils                2.9.2    2019-12-08 [2] CRAN (R 3.6.1)
    # R6                     2.4.1    2019-11-12 [2] CRAN (R 3.6.1)
    # rafalib              * 1.0.0    2015-08-09 [1] CRAN (R 3.6.1)
    # rappdirs               0.3.1    2016-03-28 [2] CRAN (R 3.6.1)
    # RColorBrewer           1.1-2    2014-12-07 [2] CRAN (R 3.6.1)
    # Rcpp                   1.0.3    2019-11-08 [2] CRAN (R 3.6.1)
    # RcppAnnoy              0.0.16   2020-03-08 [1] CRAN (R 3.6.1)
    # RCurl                  1.98-1.1 2020-01-19 [2] CRAN (R 3.6.1)
    # reshape2               1.4.3    2017-12-11 [2] CRAN (R 3.6.1)
    # rhdf5                  2.30.1   2019-11-26 [2] Bioconductor
    # Rhdf5lib               1.8.0    2019-10-29 [2] Bioconductor
    # rlang                  0.4.6    2020-05-02 [1] CRAN (R 3.6.1)
    # Rsamtools              2.2.2    2020-02-11 [2] Bioconductor
    # RSpectra               0.16-0   2019-12-01 [2] CRAN (R 3.6.1)
    # RSQLite                2.2.0    2020-01-07 [2] CRAN (R 3.6.1)
    # rsvd                   1.0.3    2020-02-17 [1] CRAN (R 3.6.1)
    # rtracklayer            1.46.0   2019-10-29 [2] Bioconductor
    # Rtsne                  0.15     2018-11-10 [1] CRAN (R 3.6.1)
    # S4Vectors            * 0.24.3   2020-01-18 [2] Bioconductor
    # scales                 1.1.0    2019-11-18 [2] CRAN (R 3.6.1)
    # scater               * 1.14.6   2019-12-16 [1] Bioconductor
    # scran                * 1.14.5   2019-11-19 [1] Bioconductor
    # segmented              1.1-0    2019-12-10 [2] CRAN (R 3.6.1)
    # sessioninfo          * 1.1.1    2018-11-05 [2] CRAN (R 3.6.1)
    # SingleCellExperiment * 1.8.0    2019-10-29 [1] Bioconductor
    # statmod                1.4.34   2020-02-17 [2] CRAN (R 3.6.1)
    # stringi                1.4.6    2020-02-17 [2] CRAN (R 3.6.1)
    # stringr                1.4.0    2019-02-10 [2] CRAN (R 3.6.1)
    # SummarizedExperiment * 1.16.1   2019-12-19 [2] Bioconductor
    # tibble                 2.1.3    2019-06-06 [2] CRAN (R 3.6.1)
    # tidyselect             1.0.0    2020-01-27 [2] CRAN (R 3.6.1)
    # uwot                   0.1.8    2020-03-16 [1] CRAN (R 3.6.1)
    # vctrs                  0.3.0    2020-05-11 [1] CRAN (R 3.6.1)
    # vipor                  0.4.5    2017-03-22 [1] CRAN (R 3.6.1)
    # viridis                0.5.1    2018-03-29 [2] CRAN (R 3.6.1)
    # viridisLite            0.3.0    2018-02-01 [2] CRAN (R 3.6.1)
    # withr                  2.1.2    2018-03-15 [2] CRAN (R 3.6.1)
    # XML                    3.99-0.3 2020-01-20 [2] CRAN (R 3.6.1)
    # XVector                0.26.0   2019-10-29 [2] Bioconductor
    # zlibbioc               1.32.0   2019-10-29 [2] Bioconductor
    # 
    # [1] /users/ntranngu/R/3.6.x
    # [2] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/site-library
    # [3] /jhpce/shared/jhpce/core/conda/miniconda3-4.6.14/envs/svnR-3.6.x/R/3.6.x/lib64/R/library

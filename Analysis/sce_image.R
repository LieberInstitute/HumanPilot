# screen -S image
# qrsh -l bluejay,mem_free=30G,h_vmem=30G,h_fsize=100G -pe local 4
# module load conda_R/3.6.x

library('SingleCellExperiment')
library('ClusterR')
library('BiocParallel')
library('ggplot2')
library('cowplot')
library('sessioninfo')

dir.create('pdf_image', showWarnings = FALSE)
dir.create('rda_image', showWarnings = FALSE)


## From convert_sce.R
load('geom_spatial.Rdata', verbose = TRUE)

## From sce_scran.R
load('Human_DLPFC_Visium_processedData_sce_scran.Rdata', verbose = TRUE)

## Using the top genes
info_full <- assays(sce)$logcounts[top.hvgs, ]
sce_df <- as.data.frame(t(as.matrix(info_full)))

## Add the row and col data, as well as the subject/position factor
sce_df$row <- unlist(lapply(split(sce$row, sce$sample_name), scale))
sce_df$col <- unlist(lapply(split(sce$col, sce$sample_name), scale))
sce_df$subject_position <- scale(as.integer(factor(sce$subject_position)))

## For making it easier later
extra_cols <- seq_len(3) + length(top.hvgs)

Sys.time()
km_init_k6_full <- KMeans_rcpp(sce_df, clusters = 6, num_init = 5, max_iters = 100,          
    initializer = 'kmeans++', verbose = TRUE, seed = 20191112)
Sys.time()

options(width = 200)
addmargins(table(sort_clusters(km_init_k6_full$clusters), sce$sample_name))

#     151507 151508 151509 151510 151669 151670 151671 151672 151673 151674 151675 151676   Sum
# 1     2138   2168   2515   2724   3011   2981   3147   3249   2567   2025   2406   2341 31272
# 2     1104   1153   1029    728    318    155    853    653    297    955    423    357  8025
# 3      876    988   1097   1037    323    349    106     88    138    108    168    169  5447
# 4      107     74    148    145      9     13      4     25    637    585    595    593  2935
# 5        0      1      0      0      0      0      0      0      0      0      0      0     1
# 6        1      0      0      0      0      0      0      0      0      0      0      0     1
# Sum   4226   4384   4789   4634   3661   3498   4110   4015   3639   3673   3592   3460 47681

addmargins(table(sort_clusters(km_init_k6_full$clusters), sce$subject_position))

#     Br5292_pos0 Br5292_pos300 Br5595_pos0 Br5595_pos300 Br8100_pos0 Br8100_pos300   Sum
# 1          4306          5239        5992          6396        4592          4747 31272
# 2          2257          1757         473          1506        1252           780  8025
# 3          1864          2134         672           194         246           337  5447
# 4           181           293          22            29        1222          1188  2935
# 5             1             0           0             0           0             0     1
# 6             1             0           0             0           0             0     1
# Sum        8610          9423        7159          8125        7312          7052 47681

addmargins(table(sort_clusters(km_init_k6_full$clusters), sce$subject))
#     Br5292 Br5595 Br8100   Sum
# 1     9545  12388   9339 31272
# 2     4014   1979   2032  8025
# 3     3998    866    583  5447
# 4      474     51   2410  2935
# 5        1      0      0     1
# 6        1      0      0     1
# Sum  18033  15284  14364 47681


sce_image_grid(sce, km_init_k6_full$clusters, 'pdf_image/grid_k6_full.pdf')
save(km_init_k6_full, file = 'rda_image/km_init_k6_full.Rdata')





Sys.time()
km_init_k6_full_noXY <- KMeans_rcpp(sce_df[, - extra_cols], clusters = 6, num_init = 5, max_iters = 100,          
    initializer = 'kmeans++', verbose = TRUE, seed = 20191112)
Sys.time()
sce_image_grid(sce, km_init_k6_full_noXY$clusters, 'pdf_image/grid_k6_full_noXY.pdf')
save(km_init_k6_full_noXY, file = 'rda_image/km_init_k6_full_noXY.Rdata')

table(
    sort_clusters(km_init_k6_full$clusters),
    sort_clusters(km_init_k6_full_noXY$clusters)
)
#       1     2     3     4     5     6
# 1 14393 11447  5426     0     1     5
# 2    16   213  5562   136    41  2057
# 3    24   558   133  4589     0   143
# 4     2     1     4     4  2875    49
# 5     0     0     0     1     0     0
# 6     0     0     0     1     0     0

addmargins(table(sort_clusters(km_init_k6_full_noXY$clusters), sce$subject))
#     Br5292 Br5595 Br8100   Sum
# 1     6198   2219   6018 14435
# 2      919   9811   1489 12219
# 3     5721   1826   3578 11125
# 4     3700    517    514  4731
# 5      445     51   2421  2917
# 6     1050    860    344  2254
# Sum  18033  15284  14364 47681


Sys.time()
km_init_k6_full_noXY_center <- KMeans_rcpp(scale(sce_df[, - extra_cols], scale = FALSE), clusters = 6, num_init = 5, max_iters = 100,          
    initializer = 'kmeans++', verbose = TRUE, seed = 20191112)
Sys.time()
sce_image_grid(sce, km_init_k6_full_noXY_center$clusters, 'pdf_image/grid_k6_full_noXY_center.pdf')
save(km_init_k6_full_noXY_center, file = 'rda_image/km_init_k6_full_noXY_center.Rdata')

table(
    sort_clusters(km_init_k6_full_noXY$clusters),
    sort_clusters(km_init_k6_full_noXY_center$clusters)
)
## centering had no effect
#       1     2     3     4     5     6
# 1 14435     0     0     0     0     0
# 2     0 12219     0     0     0     0
# 3     0     0 11125     0     0     0
# 4     0     0     0  4731     0     0
# 5     0     0     0     0  2917     0
# 6     0     0     0     0     0  2254


## No need for this one
# Sys.time()
# km_init_k6_full_center <- KMeans_rcpp(cbind(scale(sce_df[, - extra_cols], scale = FALSE), sce_df[, extra_cols]), clusters = 6, num_init = 5, max_iters = 100,
#     initializer = 'kmeans++', verbose = TRUE, seed = 20191112)
# Sys.time()
#
# sce_image_grid(sce, km_init_k6_full_center$clusters, 'pdf_image/grid_k6_full_center.pdf')



Sys.time()
km_init_k6_full_noSubjPos <- KMeans_rcpp(sce_df[, - ncol(sce_df)], clusters = 6, num_init = 5, max_iters = 100,          
    initializer = 'kmeans++', verbose = TRUE, seed = 20191112)
Sys.time()
sce_image_grid(sce, km_init_k6_full_noSubjPos$clusters, 'pdf_image/grid_k6_full_noSubjPos.pdf')
save(km_init_k6_full_noSubjPos, file = 'rda_image/km_init_k6_full_noSubjPos.Rdata')

table(
    sort_clusters(km_init_k6_full$clusters),
    sort_clusters(km_init_k6_full_noSubjPos$clusters)
)
#       1     2     3     4     5     6
# 1 13197 11903  6167     0     0     5
# 2    43    47  5687   120    23  2105
# 3     0   821   135  4363     0   128
# 4     7     0    11     2  2864    51
# 5     0     0     0     1     0     0
# 6     0     0     0     1     0     0

table(
    sort_clusters(km_init_k6_full_noXY$clusters),
    sort_clusters(km_init_k6_full_noSubjPos$clusters)
)
#       1     2     3     4     5     6
# 1 10320  4045    69     0     0     1
# 2  2117  7994  2081    23     2     2
# 3   805   476  9774    10     0    60
# 4     0   253    39  4418     2    19
# 5     5     0    26     0  2883     3
# 6     0     3    11    36     0  2204

addmargins(table(sort_clusters(km_init_k6_full_noSubjPos$clusters), sce$subject))
#     Br5292 Br5595 Br8100   Sum
# 1     3743   3669   5835 13247
# 2     4100   6579   2092 12771
# 3     5124   3655   3221 12000
# 4     3548    466    473  4487
# 5      441     50   2396  2887
# 6     1077    865    347  2289
# Sum  18033  15284  14364 47681



## Try with k-means on the 50 PCs
pcs <- reducedDim(sce, 'PCA')
dim(pcs)
# [1] 47681    50

Sys.time()
## Since it's a lot faster, I can run it with num_init = 20 instead of 5
km_init_k6_pcs_only <- KMeans_rcpp(pcs, clusters = 6, num_init = 20, max_iters = 100,          
    initializer = 'kmeans++', verbose = TRUE, seed = 20191112)
Sys.time()
## Takes about 2 minutes
# [1] "2019-11-13 15:59:58 EST"
# [1] "2019-11-13 16:01:22 EST"
sce_image_grid(sce, km_init_k6_pcs_only$clusters, 'pdf_image/grid_k6_pcs_only.pdf')
save(km_init_k6_pcs_only, file = 'rda_image/km_init_k6_pcs_only.Rdata')


addmargins(table(sort_clusters(km_init_k6_pcs_only$clusters), sce$subject))
#     Br5292 Br5595 Br8100   Sum
# 1     6277   2092   6020 14389
# 2      749   9951   1428 12128
# 3     5752   1790   3622 11164
# 4     3751    538    529  4818
# 5      448     50   2421  2919
# 6     1056    863    344  2263
# Sum  18033  15284  14364 47681

## Expression only vs PCs only: they mostly agree
table(
    sort_clusters(km_init_k6_full_noXY$clusters),
    sort_clusters(km_init_k6_pcs_only$clusters)
)
#       1     2     3     4     5     6
# 1 14150   171   111     3     0     0
# 2   202 11902    42    72     0     1
# 3    35    53 10998    19     3    17
# 4     2     1     1  4719     0     8
# 5     0     1     1     0  2914     1
# 6     0     0    11     5     2  2236


km_init_k6_pcs_noSubj <- KMeans_rcpp(cbind(pcs, sce_df[, extra_cols[-3]]), clusters = 6, num_init = 20, max_iters = 100,          
    initializer = 'kmeans++', verbose = TRUE, seed = 20191112)
sce_image_grid(sce, km_init_k6_pcs_noSubj$clusters, 'pdf_image/grid_k6_pcs_noSubj.pdf')
save(km_init_k6_pcs_noSubj, file = 'rda_image/km_init_k6_pcs_noSubj.Rdata')
addmargins(table(sort_clusters(km_init_k6_pcs_noSubj$clusters), sce$subject))
#     Br5292 Br5595 Br8100   Sum
# 1     4155   3406   5798 13359
# 2     3722   6767   2028 12517
# 3     5056   3738   3322 12116
# 4     3546    451    463  4460
# 5      444     49   2395  2888
# 6     1110    873    358  2341
# Sum  18033  15284  14364 47681

table(
    sort_clusters(km_init_k6_pcs_only$clusters),
    sort_clusters(km_init_k6_pcs_noSubj$clusters)
)
#       1     2     3     4     5     6
# 1 10629  3704    55     0     0     1
# 2  1962  7972  2180    10     2     2
# 3   765   491  9804    15     0    89
# 4     0   345    39  4407     3    24
# 5     3     0    31     0  2883     2
# 6     0     5     7    28     0  2223

## Expression vs PCs (with XY but no subj): they mostly agree
table(
    sort_clusters(km_init_k6_full_noSubjPos$clusters),
    sort_clusters(km_init_k6_pcs_noSubj$clusters)
)
#       1     2     3     4     5     6
# 1 12727   169   348     0     3     0
# 2   555 12191    18     7     0     0
# 3    77   123 11743    14     1    42
# 4     0    33     0  4438     1    15
# 5     0     1     5     0  2881     0
# 6     0     0     2     1     2  2284

km_init_k6_pcs_wXYS <- KMeans_rcpp(cbind(pcs, sce_df[, extra_cols]), clusters = 6, num_init = 20, max_iters = 100,          
    initializer = 'kmeans++', verbose = TRUE, seed = 20191112)
sce_image_grid(sce, km_init_k6_pcs_wXYS$clusters, 'pdf_image/grid_k6_pcs_wXYS.pdf')
save(km_init_k6_pcs_wXYS, file = 'rda_image/km_init_k6_pcs_wXYS.Rdata')
addmargins(table(sort_clusters(km_init_k6_pcs_wXYS$clusters), sce$subject))
#     Br5292 Br5595 Br8100   Sum
# 1     2887   3920   6278 13085
# 2     4530   6407   1904 12841
# 3     5521   3564   2957 12042
# 4     3547    459    442  4448
# 5      440     54   2414  2908
# 6     1108    880    369  2357
# Sum  18033  15284  14364 47681

table(
    sort_clusters(km_init_k6_pcs_only$clusters),
    sort_clusters(km_init_k6_pcs_wXYS$clusters)
)
#      1    2    3    4    5    6
# 1 9541 4595  252    0    0    1
# 2 2666 7426 2015   15    4    2
# 3  871  472 9717    9    1   94
# 4    0  343   46 4397    3   29
# 5    7    0    8    0 2899    5
# 6    0    5    4   27    1 2226

table(
    sort_clusters(km_init_k6_pcs_noSubj$clusters),
    sort_clusters(km_init_k6_pcs_wXYS$clusters)
)
#       1     2     3     4     5     6
# 1 11951   955   453     0     0     0
# 2   489 11858   146    22     1     1
# 3   641    12 11417     2    26    18
# 4     0    16    15  4423     0     6
# 5     4     0     1     0  2880     3
# 6     0     0    10     1     1  2329


## Expression vs PCs (with XY and subj/position info): contrary to the rest,
## these don't agree mutch
table(
    sort_clusters(km_init_k6_full$clusters),
    sort_clusters(km_init_k6_pcs_wXYS$clusters)
)
#       1     2     3     4     5     6
# 1 13030 11951  6286     0     0     5
# 2    46    49  5615   110    37  2168
# 3     0   841   137  4335     0   134
# 4     9     0     4     1  2871    50
# 5     0     0     0     1     0     0
# 6     0     0     0     1     0     0


## Explore a range of ks
ks <- 4:24

## (same colors as sce_scran.R for SNN with K = 50)
## Needs 28 colors!
## From https://medialab.github.io/iwanthue/ with the default preset
cols <- c("#de84b0",
"#78bb40",
"#9c45bd",
"#46c06f",
"#cb3e97",
"#488733",
"#b978e9",
"#beae36",
"#5c6ade",
"#db9443",
"#5985dc",
"#cf4e32",
"#43c4c4",
"#d84068",
"#5fb88e",
"#e471d1",
"#327e58",
"#7454b1",
"#a4b266",
"#964f95",
"#72722a",
"#c18cd3",
"#a06332",
"#54a4d6",
"#dc8074",
"#5465a4",
"#9f4765",
"#a09cdf")
names(cols) <- seq_len(length(cols))

## Expression only
pcs_klist_exprOnly <- lapply(ks, function(k) {
    clus_res <- KMeans_rcpp(pcs, clusters = k,
        num_init = 10, max_iters = 100,          
        initializer = 'kmeans++', verbose = FALSE, seed = 20191112)
    sort_clusters(clus_res$clusters)
})

pdf('pdf_image/grid_kmeans_pcs_exprOnly.pdf', height = 24, width = 36)
mapply(function(k, clus) {
    plots <- sce_image_grid(sce, clus,
        ... = paste(' k = ', k, 'expr PCs only'),
        return_plots = TRUE, colors = cols)
    print(cowplot::plot_grid(plotlist = plots))
    return(invisible(NULL))
}, ks, pcs_klist_exprOnly)
dev.off()
save(pcs_klist_exprOnly, file = 'rda_image/pcs_klist_exprOnly.Rdata')


## with X and Y
pcs_klist_withXY <- lapply(ks, function(k) {
    clus_res <- KMeans_rcpp(cbind(pcs, sce_df[, extra_cols[-3]]), clusters = k,
        num_init = 10, max_iters = 100,          
        initializer = 'kmeans++', verbose = FALSE, seed = 20191112)
    sort_clusters(clus_res$clusters)
})

pdf('pdf_image/grid_kmeans_pcs_withXY.pdf', height = 24, width = 36)
mapply(function(k, clus) {
    plots <- sce_image_grid(sce, clus,
        ... = paste(' k = ', k, 'expr PCs with XY'),
        return_plots = TRUE, colors = cols)
    print(cowplot::plot_grid(plotlist = plots))
    return(invisible(NULL))
}, ks, pcs_klist_withXY)
dev.off()
save(pcs_klist_withXY, file = 'rda_image/pcs_klist_withXY.Rdata')


## with X and Y & subj/position
pcs_klist_withXY_andSubj <- lapply(ks, function(k) {
    clus_res <- KMeans_rcpp(cbind(pcs, sce_df[, extra_cols]), clusters = k,
        num_init = 10, max_iters = 100,          
        initializer = 'kmeans++', verbose = FALSE, seed = 20191112)
    sort_clusters(clus_res$clusters)
})

pdf('pdf_image/grid_kmeans_pcs_withXY_andSubj.pdf', height = 24, width = 36)
mapply(function(k, clus) {
    plots <- sce_image_grid(sce, clus,
        ... = paste(' k = ', k, 'expr PCs with XY & subj/pos'),
        return_plots = TRUE, colors = cols)
    print(cowplot::plot_grid(plotlist = plots))
    return(invisible(NULL))
}, ks, pcs_klist_withXY_andSubj)
dev.off()
save(pcs_klist_withXY_andSubj, file = 'rda_image/pcs_klist_withXY_andSubj.Rdata')


## Try a wild idea
pcs_klist_exprOnly_mat <- scale(do.call(cbind, pcs_klist_exprOnly))

pcs_klist_exprOnly_on_ks <- lapply(3:(length(pcs_klist_exprOnly)-1), function(k) {
    clus_res <- KMeans_rcpp(pcs_klist_exprOnly_mat, clusters = k,
        num_init = 10, max_iters = 100,          
        initializer = 'kmeans++', verbose = FALSE, seed = 20191112)
    sort_clusters(clus_res$clusters)
})

pdf('pdf_image/grid_kmeans_pcs_exprOnly_on_ks.pdf', height = 24, width = 36)
mapply(function(k, clus) {
    plots <- sce_image_grid(sce, clus,
        ... = paste(' k = ', k, "on k's expr PCs only"),
        return_plots = TRUE, colors = cols)
    print(cowplot::plot_grid(plotlist = plots))
    return(invisible(NULL))
}, 3:(length(pcs_klist_exprOnly)-1), pcs_klist_exprOnly_on_ks)
dev.off()
save(pcs_klist_exprOnly_on_ks, file = 'rda_image/pc_klist_exprOnly_on_ks.Rdata')


pcs_klist_exprOnly_on_ks_withXY <- lapply(3:(length(pcs_klist_exprOnly)), function(k) {
    clus_res <- KMeans_rcpp(cbind(pcs_klist_exprOnly_mat, sce_df[, extra_cols[-3]]), clusters = k,
        num_init = 10, max_iters = 100,          
        initializer = 'kmeans++', verbose = FALSE, seed = 20191112)
    sort_clusters(clus_res$clusters)
})

pdf('pdf_image/grid_kmeans_pcs_exprOnly_on_ks_withXY.pdf', height = 24, width = 36)
mapply(function(k, clus) {
    plots <- sce_image_grid(sce, clus,
        ... = paste(' k = ', k, "on k's & XY on expr PCs only"),
        return_plots = TRUE, colors = cols)
    print(cowplot::plot_grid(plotlist = plots))
    return(invisible(NULL))
}, 3:(length(pcs_klist_exprOnly)), pcs_klist_exprOnly_on_ks_withXY)
dev.off()
save(pcs_klist_exprOnly_on_ks_withXY, file = 'rda_image/pc_klist_exprOnly_on_ks_withXY.Rdata')


pcs_klist_exprOnly_on_ks_withXY_andSubj <- lapply(3:(length(pcs_klist_exprOnly)+1), function(k) {
    clus_res <- KMeans_rcpp(cbind(pcs_klist_exprOnly_mat, sce_df[, extra_cols]), clusters = k,
        num_init = 10, max_iters = 100,          
        initializer = 'kmeans++', verbose = FALSE, seed = 20191112)
    sort_clusters(clus_res$clusters)
})

pdf('pdf_image/grid_kmeans_pcs_exprOnly_on_ks_withXY_andSubj.pdf', height = 24, width = 36)
mapply(function(k, clus) {
    plots <- sce_image_grid(sce, clus,
        ... = paste(' k = ', k, "on k's & XY + subj/pos on expr PCs only"),
        return_plots = TRUE, colors = cols)
    print(cowplot::plot_grid(plotlist = plots))
    return(invisible(NULL))
}, 3:(length(pcs_klist_exprOnly)+1), pcs_klist_exprOnly_on_ks_withXY_andSubj)
dev.off()
save(pcs_klist_exprOnly_on_ks_withXY_andSubj, file = 'rda_image/pc_klist_exprOnly_on_ks_withXY_andSubj.Rdata')



pdf('pdf_image/grid_kmeans_pcs_exprOnly_VS_pcs_exprOnly_on_ks.pdf', height = 24, width = 36)
mapply(function(k, clus1, clus2) {
    plots <- sce_image_grid_comp(sce, clus1, clus2,
        map_subset = sce$subject == 'Br5292',
        ... = paste(' k = ', k, "on k's & XY + subj/pos on expr PCs only vs expr PCs only"),
        return_plots = TRUE)
    print(cowplot::plot_grid(plotlist = plots))
    return(invisible(NULL))
}, 4:8, pcs_klist_exprOnly[1:5], pcs_klist_exprOnly_on_ks[-1])
dev.off()



clus1 <- pcs_klist_exprOnly[[2]]
clus2 <- pcs_klist_exprOnly_on_ks[[3]]
addmargins(table(clus1, clus2))

addmargins(table(sort_clusters(clus1, map_subset = sce$subject == 'Br5292'), sort_clusters(clus2, map_subset = sce$subject == 'Br5292'))) 
 

### Stuff that didn't work


# sce <- sce[top.hvgs, sce$sample_name == '151673']
# slotname <- 'logcounts'
# ## To array format
# sce_to_array <- function(sce, slotname) {
#     info <- assays(sce)[[slotname]]
#
#     rows <- sort(unique(sce$row))
#     cols <- sort(unique(sce$col))
#     ids <- paste0(
#         rep(rows, each = length(cols)), '_',
#         rep(cols, length(rows))
#     )
#
#     ids_obs <- paste0(sce$row, '_', sce$col)
#     m <- match(ids, ids_obs)
#     table(is.na(m))
#
#     arr <- array(NA, c(length(rows), length(cols), nrow(sce)))
#
#     ids_pairs <- ids[match(ids_obs, ids)]
#     ids_pairs_i <- match(gsub('_.*', '', ids_pairs), rows)
#     ids_pairs_j <- match(gsub('.*_', '', ids_pairs), cols)
#
#     for(g in seq_len(ncol(sce))) {
#         arr[ids_pairs_i[g], ids_pairs_j[g], ] <- info[, g]
#     }
#
#     arr_vec <- apply(arr, 3, as.vector)
#     arr_vec <- arr_vec[complete.cases(arr_vec), ]
#
#     return(arr_vec)
# }

## Failed completely
# km_mb <- MiniBatchKmeans(arr_vec, clusters = 9, batch_size = 20, num_init = 5, max_iters = 100,
#     init_fraction = 0.2, initializer = 'kmeans++', early_stop_iter = 10,
#     verbose = FALSE)
# pr_mb <- predict_MBatchKMeans(arr_vec, km_mb$centroids)
# sce$km_mb_cluster <- as.vector(pr_mb)[m[!is.na(m)]]
# pdf('pdf_image/test_151673_k9.pdf')
# sce_image_clus(sce, '151673', 'km_mb_cluster')
# dev.off()


# library('tidyr')
#
# sce_to_df <- function(sce, slotname) {
#     info <- assays(sce)[[slotname]]
#
#     info <- as.data.frame(t(as.matrix(info)))
#     colnames(info) <- seq_len(nrow(sce))
#
#     info$row <- sce$row
#     info$col <- sce$col
#     sce_df <- gather(info, key = 'gene', value = 'expr', -(row:col))
#     sce_df$umi <- rep(seq_len(ncol(sce)), nrow(sce))
#     sce_df$gene <- as.integer(sce_df$gene)
#
#     return(sce_df)
# }


# sce_df_scaled <- scale(sce_df)
# summary(sce_df_scaled)
#
# set.seed(20191112)
# ## Takes about 2-3 min
# km_mb_df <- MiniBatchKmeans(sce_df_scaled, clusters = 9, batch_size = 20,
#     num_init = 5, max_iters = 100,
#     init_fraction = 0.2, initializer = 'kmeans++', early_stop_iter = 10,
#     verbose = FALSE)
# pr_mb_df <- predict_MBatchKMeans(sce_df_scaled, km_mb_df$centroids)
#
# table(pr_mb_df)
# #      1       2       3       4       5       6       7       8       9
# # 979814   31770  728867  576061  929346  602235 1686166  856007  808358
#
# table(pr_mb_df[seq_len(ncol(sce))])
# #  1    2    3    4    5    6    7    9
# # 51  309  406 1633    1  950  110  172
#
# cl_umi <- table(sce_df$umi, pr_mb_df)
#
# sce$clus_k9_df <- apply(cl_umi, 1, which.max)
# table(sce$clus_k9_df)
# #   1    3    5    6    7    8    9
# # 595  416  556    8 1113  609  335
#
# pdf('pdf_image/test_151673_k9_df.pdf')
# sce_image_clus(sce, '151673', 'clus_k9_df')
# dev.off()


## This seemed to work (only using sample 9)

sce_to_df <- function(sce, slotname) {
    info <- assays(sce)[[slotname]]

    sce_df <- as.data.frame(t(as.matrix(info)))

    sce_df$row <- sce$row
    sce_df$col <- sce$col

    return(sce_df)
}


sce_df_scaled <- scale(sce_df)

# d <- dist(sce_df_scaled)
# km <- kmeans(d, centers = 9, iter.max = 100)

km_init <- KMeans_rcpp(sce_df_scaled, clusters = 9, num_init = 20, max_iters = 100,          
    initializer = 'kmeans++', verbose = TRUE, seed = 20191112)

table(km_init$clusters)
# 1    2    3    4    5    6    7    8    9
# 1    1  287    4    4 2217  754  363    1

sce$clus_k9_df <- km_init$clusters
pdf('pdf_image/test_151673_k9_df_take2.pdf')
sce_image_clus(sce, '151673', 'clus_k9_df')
dev.off()



km_init_k6 <- KMeans_rcpp(sce_df_scaled, clusters = 6, num_init = 20, max_iters = 100,          
    initializer = 'kmeans++', verbose = TRUE, seed = 20191112)

table(km_init_k6$clusters)
# 1    2    3    4    5    6
# 1  314  360    1  935 2021

sce$clus_k6_df <- km_init_k6$clusters
pdf('pdf_image/test_151673_k6_df_take2.pdf')
sce_image_clus(sce, '151673', 'clus_k6_df')
dev.off()


km_init_k6_noXY <- KMeans_rcpp(sce_df_scaled[, -c(1983, 1984)], clusters = 6, num_init = 20, max_iters = 100,          
    initializer = 'kmeans++', verbose = TRUE, seed = 20191112)

table(km_init_k6_noXY$clusters)
# 1    2    3    4    5    6
# 1  309  364    1  915 2042

table(km_init_k6_noXY$clusters, km_init_k6$clusters)
#      1    2    3    4    5    6
# 1    1    0    0    0    0    0
# 2    0  309    0    0    0    0
# 3    0    5  359    0    0    0
# 4    0    0    0    1    0    0
# 5    0    0    1    0  910    4
# 6    0    0    0    0   25 2017

sce$clus_k6_df_noXY <- km_init_k6_noXY$clusters
pdf('pdf_image/test_151673_k6_df_take2_noXY.pdf')
sce_image_clus(sce, '151673', 'clus_k6_df_noXY')
dev.off()


km_init_k6_noXY_noScale <- KMeans_rcpp(sce_df[, -c(1983, 1984)], clusters = 6, num_init = 20, max_iters = 100,          
    initializer = 'kmeans++', verbose = TRUE, seed = 20191112)

table(km_init_k6_noXY_noScale$clusters)
#   1    2    3    4    5    6
# 321    1  341  390  923 1656

table(km_init_k6_noXY_noScale$clusters, km_init_k6$clusters)
#      1    2    3    4    5    6
# 1    0    0  284    1   35    1
# 2    0    0    1    0    0    0
# 3    0    0    0    0  339    2
# 4    1  314   75    0    0    0
# 5    0    0    0    0  522  401
# 6    0    0    0    0   39 1617

table(km_init_k6_noXY_noScale$clusters, km_init_k6_noXY$clusters)
#      1    2    3    4    5    6
# 1    0    0  283    1   36    1
# 2    0    0    1    0    0    0
# 3    0    0    0    0  337    4
# 4    1  309   80    0    0    0
# 5    0    0    0    0  512  411
# 6    0    0    0    0   30 1626

sce$clus_k6_df_noXY_noScale <- km_init_k6_noXY_noScale$clusters
pdf('pdf_image/test_151673_k6_df_take2_noXY_noScale.pdf')
sce_image_clus(sce, '151673', 'clus_k6_df_noXY_noScale')
dev.off()


km_init_k6_noScale <- KMeans_rcpp(sce_df, clusters = 6, num_init = 20, max_iters = 100,          
    initializer = 'kmeans++', verbose = TRUE, seed = 20191112)

table(km_init_k6_noScale$clusters)
#   1   2   3   4   5   6
# 579 628 579 603 660 583

table(km_init_k6_noXY_noScale$clusters, km_init_k6_noScale$clusters)
#     1   2   3   4   5   6
# 1   1 151 164   2   3   0
# 2   0   0   1   0   0   0
# 3  71  13   3  48 197   9
# 4   0  14 376   0   0   0
# 5  73 325  32  57 193 243
# 6 434 125   3 496 267 331

sce$clus_k6_df_noScale <- km_init_k6_noScale$clusters
pdf('pdf_image/test_151673_k6_df_take2_noScale.pdf')
sce_image_clus(sce, '151673', 'clus_k6_df_noScale')
dev.off()


km_init_k6_noScale_sXY <- KMeans_rcpp(cbind(sce_df[, -c(1983, 1984)], sce_df_scaled[, c(1983, 1984)]), clusters = 6, num_init = 20, max_iters = 100,          
    initializer = 'kmeans++', verbose = TRUE, seed = 20191112)

table(km_init_k6_noScale_sXY$clusters)
#   1    2    3    4    5    6
# 857    1  391 1693  370  320

table(km_init_k6_noXY_noScale$clusters, km_init_k6_noScale_sXY$clusters)
#      1    2    3    4    5    6
# 1    1    0    0    0    0  320
# 2    0    0    1    0    0    0
# 3   10    0    0    0  331    0
# 4    0    0  390    0    0    0
# 5  811    1    0   87   24    0
# 6   35    0    0 1606   15    0

table(km_init_k6$clusters, km_init_k6_noScale_sXY$clusters)
#      1    2    3    4    5    6
# 1    0    0    1    0    0    0
# 2    0    0  314    0    0    0
# 3    0    0   76    0    0  284
# 4    0    0    0    0    0    1
# 5  477    1    0   67  356   34
# 6  380    0    0 1626   14    1

sce$clus_k6_df_noScale_sXY <- km_init_k6_noScale_sXY$clusters
pdf('pdf_image/test_151673_k6_df_take2_noScale_sXY.pdf')
sce_image_clus(sce, '151673', 'clus_k6_df_noScale_sXY')
dev.off()





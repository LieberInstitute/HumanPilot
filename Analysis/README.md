Analysis
========


## Code by Leo

### Summary

* [`convert_sce.R`](convert_sce.R): builds the SingleCellExperiment Bioconductor object and defines some plotting functions used in other scripts.
* [`sce_scran.R`](sce_scran.R): determines the 1942 variable genes, computes PCs and saves the `sce` object (file `Human_DLPFC_Visium_processedData_sce_scran.Rdata`) used by other scripts. Also explores several things you can do based on https://osca.bioconductor.org/ and the `scran` vignette.
* [`sce_zinbwave.R`](sce_zinbwave.R): uses the `sce` object from `sce_scran.R` and runs `zinbwave`.
* [`sce_image.R`](sce_image.R): uses the `sce` object from `sce_scran.R` and runs several k-means approaches using the `clusteR` package.

### [`convert_sce.R`](convert_sce.R)

Related files:

* `/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Human_DLPFC_Visium_processedData_sce.Rdata`
* `/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/geom_spatial.Rdata`

Details:

* Builds the `sce` object with image data under `metadata(sce)$image` which is a single data.frame. Subsetting doesn't automatically subset the image, so you have to do it yourself when plotting.
* Creates some of the functions and the ggplot2 `geom_spatial()` layer used for plotting later on. Newer and more complex versions are in the [`global.R` file from `spatialLIBD`](https://github.com/LieberInstitute/spatialLIBD/blob/master/global.R). In particular, the colors are different between these versions as now `spatialLIBD` uses the original colors from the 10X scripts when the number of clusters <= 12 and then relies on `Polychrome::palette36.colors()` (which would break if more than 36 colors are needed).

### [`sce_scran.R`](sce_scran.R)

Related files:

* [pdf_scran](pdf_scran/)
* [rda_scran](rda_scran/)
* `/dcs04/lieber/lcolladotor/with10x_LIBD001/HumanPilot/Analysis/Human_DLPFC_Visium_processedData_sce_scran.Rdata`

Details:

* Defines the `sce` and `top.hvgs` objects used in other scripts and in `spatialLIBD`. As we decide to add more stuff to the `sce` object, I'll edit this R code as well.
* Builds a SNN graph with 50 nearest neighbors (K=50) that results in 28 clusters across all 12 images.
* Cuts that SNN graph to k = 4 up to 28.
* Mostly follows https://osca.bioconductor.org/.

It's long enough now that I'll continue new analyses in other scripts.


### [`sce_zinbwave.R`](sce_zinbwave.R)

Related files:

* [pdf_zinbwave](pdf_zinbwave/)
* [rda_zinbwave](rda_zinbwave/)

Details:

* Runs `zinbwave::zinbwave()` and saves the results. Initially, it chose the genes using `zinbwave` but now it uses the same ones from `sce_scran.R`. It failed at running `RSEC()` either due to memory (when using `kmeans`) or because I stopped trying to get `ClusteR::KMeans_rcpp` to work with `RSEC()`.


### [`sce_image.R`](sce_image.R)

Related files:

* [pdf_image](pdf_image/)
* [rda_image](rda_image/)

Details:

* The name of the script comes from me thinking that we could process it like an image as in https://cran.rstudio.com/web/packages/ClusterR/vignettes/the_clusterR_package.html (the dog example). But then I realized that the image example didn't use X or Y information. Some of this quick exploration code is near the end of the script (just so I wouldn't lose the code forever).
* Ultimately, the likely useful output of this script is the set of k-means clusters from k = 4 up to 28 using `ClusteR`. Although I didn't save all the kmeans output (only the cluster labels), so we might need to re-compute it.
* This is also the script where I played the most with trying to add the X and Y information as well as a blocking factor (6 levels, one per each subject + slice (position) pair of images).


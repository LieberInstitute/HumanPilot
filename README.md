HumanPilot
==========

Code for the 10x Spatial Transcriptomics `HumanPilot` project.

# `Analysis` directory

* `convert_sce.R`: builds the SingleCellExperiment Bioconductor object and defines some plotting functions used in other scripts.
* `sce_scran.R`: determines the 1942 variable genes, computes PCs and saves the `sce` object (file `Human_DLPFC_Visium_processedData_sce_scran.Rdata`) used by other scripts. Also explores several things you can do based on https://osca.bioconductor.org/ and the `scran` vignette.
* `sce_zinbwave.R`: uses the `sce` object from `sce_scran.R` and runs `zinbwave`.
* `sce_image.R`: uses the `sce` object from `sce_scran.R` and runs several k-means approaches using the `clusteR` package.

# Internal

* JHPCE location: `/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot`
* Main `sce` R object file: `/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Human_DLPFC_Visium_processedData_sce_scran.Rdata`.

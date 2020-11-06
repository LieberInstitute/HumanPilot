
<!-- README.md is generated from README.Rmd. Please edit that file -->

# HumanPilot <img src="http://research.libd.org/spatialLIBD/reference/figures/logo.png" align="right" />

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/225910046.svg)](https://zenodo.org/badge/latestdoi/225910046)
<!-- badges: end -->

Welcome to the `spatialLIBD` project\! It is composed of the
`HumanPilot` described here as well as:

  - a [shiny](https://shiny.rstudio.com/) web application that we are
    hosting at
    [spatial.libd.org/spatialLIBD/](http://spatial.libd.org/spatialLIBD/)
    that can handle a
    [limited](https://github.com/LieberInstitute/spatialLIBD/issues/2)
    set of concurrent users,
  - a Bioconductor package at
    [bioconductor.org/packages/spatialLIBD](http://bioconductor.org/packages/spatialLIBD)
    (or from [here](http://research.libd.org/spatialLIBD/)) that lets
    you analyze the data and run a local version of our web application
    (with our data or yours),
  - and a [research
    article](https://www.biorxiv.org/content/10.1101/2020.02.28.969931v1)
    with the scientific knowledge we drew from this dataset. The
    analysis code for our project is available
    [here](https://github.com/LieberInstitute/HumanPilot/) that you are
    looking at right now.

This web application allows you to browse the LIBD human dorsolateral
pre-frontal cortex (DLPFC) spatial transcriptomics data generated with
the 10x Genomics Visium platform. Through the [R/Bioconductor
package](https://bioconductor.org/packages/spatialLIBD) you can also
download the data as well as visualize your own datasets using this web
application. Please check the [bioRxiv
pre-print](https://www.biorxiv.org/content/10.1101/2020.02.28.969931v1)
for more details about this project.

If you tweet about this website, the data or the R package please use
the <code>\#spatialLIBD</code> hashtag. You can find previous tweets
that way as shown
<a href="https://twitter.com/search?q=%23spatialLIBD&src=typed_query">here</a>.
Thank you\!

## Study design

As a quick overview, the data presented here is from portion of the
DLPFC that spans six neuronal layers plus white matter (**A**) for a
total of three subjects with two pairs of spatially adjacent replicates
(**B**). Each dissection of DLPFC was designed to span all six layers
plus white matter (**C**). Using this web application you can explore
the expression of known genes such as *SNAP25* (**D**, a neuronal gene),
*MOBP* (**E**, an oligodendrocyte gene), and known layer markers from
mouse studies such as *PCP4* (**F**, a known layer 5 marker gene).

<img src="http://research.libd.org/spatialLIBD/reference/figures/paper_figure1.jpg" align="center" width="800px" />

This web application was built such that we could annotate the spots to
layers as you can see under the **spot-level data** tab. Once we
annotated each spot to a layer, we compressed the information by a
pseudo-bulking approach into **layer-level data**. We then analyzed the
expression through a set of models whose results you can also explore
through this web application. Finally, you can upload your own gene sets
of interest as well as layer enrichment statistics and compare them with
our LIBD Human DLPFC Visium dataset.

If you are interested in running this web application locally, you can
do so thanks to the `spatialLIBD` R/Bioconductor package that powers
this web application as shown below.

``` r
## Run this web application locally
spatialLIBD::run_app()

## You will have more control about the length of the
## session and memory usage.

## You could also use this function to visualize your
## own data given some requirements described
## in detail in the package vignette documentation
## at http://research.libd.org/spatialLIBD/.
```

## Shiny website mirrors

  - [Main shiny application
    website](http://spatial.libd.org/spatialLIBD)
  - [Shinyapps](https://jhubiostatistics.shinyapps.io/spatialLIBD/)
  - [Shinyapps Mirror
    1](https://jhubiostatistics.shinyapps.io/spatialLIBD_mirror01/)
  - [Shinyapps Mirror
    2](https://jhubiostatistics.shinyapps.io/spatialLIBD_mirror02/)

## R/Bioconductor package

The `spatialLIBD` package contains functions for:

  - Accessing the spatial transcriptomics data from the LIBD Human Pilot
    project ([code on
    GitHub](https://github.com/LieberInstitute/HumanPilot)) generated
    with the Visium platform from 10x Genomics. The data is retrieved
    from [Bioconductor](http://bioconductor.org/)’s `ExperimentHub`.
  - Visualizing the spot-level spatial gene expression data and
    clusters.
  - Inspecting the data interactively either on your computer or through
    [spatial.libd.org/spatialLIBD/](http://spatial.libd.org/spatialLIBD/).

For more details, please check the [documentation
website](http://lieberinstitute.github.io/spatialLIBD) or the
Bioconductor package landing page
[here](https://bioconductor.org/packages/spatialLIBD).

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `spatialLIBD` using
from [Bioconductor](http://bioconductor.org/) the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("spatialLIBD")
```

## Access the data

Through the `spatialLIBD` package you can access the processed data in
it’s final R format. However, we also provide a table of links so you
can download the raw data we received from 10x Genomics.

### Processed data

Using `spatialLIBD` you can access the Human DLPFC spatial
transcriptomics data from the 10x Genomics Visium platform. For example,
this is the code you can use to access the layer-level data. For more
details, check the help file for `fetch_data()`.

``` r
## Load the package
library('spatialLIBD')

## Download the spot-level data
sce <- fetch_data(type = 'sce')
#> Loading objects:
#>   sce

## This is a SingleCellExperiment object
sce
#> class: SingleCellExperiment 
#> dim: 33538 47681 
#> metadata(1): image
#> assays(2): counts logcounts
#> rownames(33538): ENSG00000243485 ENSG00000237613 ... ENSG00000277475
#>   ENSG00000268674
#> rowData names(9): source type ... gene_search is_top_hvg
#> colnames(47681): AAACAACGAATAGTTC-1 AAACAAGTATCTCCCA-1 ...
#>   TTGTTTCCATACAACT-1 TTGTTTGTGTAAATTC-1
#> colData names(73): barcode sample_name ... pseudobulk_UMAP_spatial
#>   markers_UMAP_spatial
#> reducedDimNames(6): PCA TSNE_perplexity50 ... TSNE_perplexity80
#>   UMAP_neighbors15
#> spikeNames(0):
#> altExpNames(0):

## Note the memory size
pryr::object_size(sce)
#> 2.08 GB

## Remake the logo image with histology information
sce_image_clus(
    sce = sce,
    clustervar = 'layer_guess_reordered',
    sampleid = '151673',
    colors = libd_layer_colors,
    ... = ' DLPFC Human Brain Layers\nMade with github.com/LieberInstitute/spatialLIBD'
)
```

<img src="http://research.libd.org/spatialLIBD/reference/figures/README-access_data-1.png" width="600px" align="center" />

### Raw data

You can access all the raw data through
[Globus](http://research.libd.org/globus/) (`jhpce#HumanPilot10x`).
Furthermore, below you can find the links to the raw data we received
from 10x Genomics.

| SampleID | h5\_filtered                                                                                    | h5\_raw                                                                                    | image\_full                                                                          | image\_hi                                                                                    | image\_lo                                                                                     | loupe                                                                       | HTML\_report                                                                                           |
| -------: | :---------------------------------------------------------------------------------------------- | :----------------------------------------------------------------------------------------- | :----------------------------------------------------------------------------------- | :------------------------------------------------------------------------------------------- | :-------------------------------------------------------------------------------------------- | :-------------------------------------------------------------------------- | :----------------------------------------------------------------------------------------------------- |
|   151507 | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151507_filtered_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151507_raw_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151507_full_image.tif) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151507_tissue_hires_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151507_tissue_lowres_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/loupe/151507.cloupe) | [GitHub](https://github.com/LieberInstitute/HumanPilot/blob/master/10X/151507/151507_web_summary.html) |
|   151508 | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151508_filtered_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151508_raw_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151508_full_image.tif) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151508_tissue_hires_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151508_tissue_lowres_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/loupe/151508.cloupe) | [GitHub](https://github.com/LieberInstitute/HumanPilot/blob/master/10X/151508/151508_web_summary.html) |
|   151509 | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151509_filtered_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151509_raw_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151509_full_image.tif) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151509_tissue_hires_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151509_tissue_lowres_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/loupe/151509.cloupe) | [GitHub](https://github.com/LieberInstitute/HumanPilot/blob/master/10X/151509/151509_web_summary.html) |
|   151510 | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151510_filtered_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151510_raw_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151510_full_image.tif) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151510_tissue_hires_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151510_tissue_lowres_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/loupe/151510.cloupe) | [GitHub](https://github.com/LieberInstitute/HumanPilot/blob/master/10X/151510/151510_web_summary.html) |
|   151669 | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151669_filtered_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151669_raw_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151669_full_image.tif) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151669_tissue_hires_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151669_tissue_lowres_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/loupe/151669.cloupe) | [GitHub](https://github.com/LieberInstitute/HumanPilot/blob/master/10X/151669/151669_web_summary.html) |
|   151670 | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151670_filtered_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151670_raw_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151670_full_image.tif) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151670_tissue_hires_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151670_tissue_lowres_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/loupe/151670.cloupe) | [GitHub](https://github.com/LieberInstitute/HumanPilot/blob/master/10X/151670/151670_web_summary.html) |
|   151671 | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151671_filtered_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151671_raw_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151671_full_image.tif) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151671_tissue_hires_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151671_tissue_lowres_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/loupe/151671.cloupe) | [GitHub](https://github.com/LieberInstitute/HumanPilot/blob/master/10X/151671/151671_web_summary.html) |
|   151672 | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151672_filtered_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151672_raw_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151672_full_image.tif) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151672_tissue_hires_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151672_tissue_lowres_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/loupe/151672.cloupe) | [GitHub](https://github.com/LieberInstitute/HumanPilot/blob/master/10X/151672/151672_web_summary.html) |
|   151673 | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151673_filtered_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151673_raw_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151673_full_image.tif) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151673_tissue_hires_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151673_tissue_lowres_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/loupe/151673.cloupe) | [GitHub](https://github.com/LieberInstitute/HumanPilot/blob/master/10X/151673/151673_web_summary.html) |
|   151674 | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151674_filtered_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151674_raw_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151674_full_image.tif) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151674_tissue_hires_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151674_tissue_lowres_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/loupe/151674.cloupe) | [GitHub](https://github.com/LieberInstitute/HumanPilot/blob/master/10X/151674/151674_web_summary.html) |
|   151675 | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151675_filtered_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151675_raw_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151675_full_image.tif) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151675_tissue_hires_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151675_tissue_lowres_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/loupe/151675.cloupe) | [GitHub](https://github.com/LieberInstitute/HumanPilot/blob/master/10X/151675/151675_web_summary.html) |
|   151676 | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151676_filtered_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151676_raw_feature_bc_matrix.h5) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151676_full_image.tif) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151676_tissue_hires_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/images/151676_tissue_lowres_image.png) | [AWS](https://spatial-dlpfc.s3.us-east-2.amazonaws.com/loupe/151676.cloupe) | [GitHub](https://github.com/LieberInstitute/HumanPilot/blob/master/10X/151676/151676_web_summary.html) |

## Citation

Below is the citation output from using `citation('spatialLIBD')` in R.
Please run this yourself to check for any updates on how to cite
**spatialLIBD**.

``` r
citation('spatialLIBD')
#> 
#> Collado-Torres L, Maynard KR, Jaffe AE (2020). _LIBD Visium spatial
#> transcriptomics human pilot data inspector_. doi:
#> 10.18129/B9.bioc.spatialLIBD (URL:
#> https://doi.org/10.18129/B9.bioc.spatialLIBD),
#> https://github.com/LieberInstitute/spatialLIBD - R package version
#> 1.0.0, <URL: http://www.bioconductor.org/packages/spatialLIBD>.
#> 
#> Maynard KR, Collado-Torres L, Weber LM, Uytingco C, Barry BK, Williams
#> SR, II JLC, Tran MN, Besich Z, Tippani M, Chew J, Yin Y, Kleinman JE,
#> Hyde TM, Rao N, Hicks SC, Martinowich K, Jaffe AE (2020).
#> "Transcriptome-scale spatial gene expression in the human dorsolateral
#> prefrontal cortex." _bioRxiv_. doi: 10.1101/2020.02.28.969931 (URL:
#> https://doi.org/10.1101/2020.02.28.969931), <URL:
#> https://www.biorxiv.org/content/10.1101/2020.02.28.969931v1>.
#> 
#> To see these entries in BibTeX format, use 'print(<citation>,
#> bibtex=TRUE)', 'toBibtex(.)', or set
#> 'options(citation.bibtex.max=999)'.
```

# `HumanPilot` code

## Re-shaping your data to our structure

As described in the `spatialLIBD` vignette, you can see the scripts in
this repository for re-shaping your data to look like ours. That is.

  - `reorganize_folder.R` available
    [here](https://github.com/LieberInstitute/HumanPilot/blob/master/reorganize_folder.R)
    re-organizes the raw data we were sent by 10x Genomics.
  - `Layer_Notebook.R` available
    [here](https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/Layer_Notebook.R)
    reads in the Visium data and builds a list of
    `RangeSummarizedExperiment()` objects from
    *[SummarizedExperiment](https://bioconductor.org/packages/3.11/SummarizedExperiment)*,
    one per sample (image) that is eventually saved as
    `Human_DLPFC_Visium_processedData_rseList.rda`.
  - `convert_sce.R` available
    [here](https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/convert_sce.R)
    reads in `Human_DLPFC_Visium_processedData_rseList.rda` and builds
    an initial `sce` object with image data under `metadata(sce)$image`
    which is a single data.frame. Subsetting doesn’t automatically
    subset the image, so you have to do it yourself when plotting as is
    done by `sce_image_clus_p()` and `sce_image_gene_p()`. Having the
    data from all images in a single object allows you to use the
    spot-level data from all images to compute clusters and do other
    similar analyses to the ones you would do with sc/snRNA-seq data.
    The script creates the `Human_DLPFC_Visium_processedData_sce.Rdata`
    file.
  - `sce_scran.R` available
    [here](https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/sce_scran.R)
    then uses *[scran](https://bioconductor.org/packages/3.11/scran)* to
    read in `Human_DLPFC_Visium_processedData_sce.Rdata`, compute the
    highly variable genes (stored in our final `sce` object at
    `rowData(sce)$is_top_hvg`), perform dimensionality reduction (PCA,
    TSNE, UMAP) and identify clusters using the data from all images.
    The resulting data is then stored as
    `Human_DLPFC_Visium_processedData_sce_scran.Rdata` and is the main
    object used throughout our analysis code
    <a id='cite-Maynard_2020'></a>(<a href='https://www.biorxiv.org/content/10.1101/2020.02.28.969931v1'>Maynard,
    Collado-Torres, Weber, Uytingco, et al., 2020</a>).
  - `make-data_spatialLIBD.R` available in the source version of
    `spatialLIBD` and [online
    here](https://github.com/LieberInstitute/spatialLIBD/blob/master/inst/scripts/make-data_spatialLIBD.R)
    is the script that reads in
    `Human_DLPFC_Visium_processedData_sce_scran.Rdata` as well as some
    other outputs from our analysis and combines them into the final
    `sce` and `sce_layer` objects provided by
    *[spatialLIBD](https://bioconductor.org/packages/3.11/spatialLIBD)*
    <a id='cite-Collado-Torres_2020'></a>(<a href='http://www.bioconductor.org/packages/spatialLIBD'>Collado-Torres,
    Maynard, and Jaffe, 2020</a>). This script simplifies some
    operations in order to simplify the code behind the
    *[shiny](https://CRAN.R-project.org/package=shiny)* application
    provided by
    *[spatialLIBD](https://bioconductor.org/packages/3.11/spatialLIBD)*.

## [10X](10X/) directory

Contains some of the raw files provided by 10X. Given their size, we
only included the small ones here. \#\# [Analysis](Analysis/) directory

The `README.md` was the one we initially prepared for our collaborators
at an early stage of the project. That README file described some of our
initial explorations using packages such as
*[scran](https://bioconductor.org/packages/3.11/scran)*,
*[zinbwave](https://bioconductor.org/packages/3.11/zinbwave)* and other
approaches such as using k-means with X/Y spatial information. These
analyses were not used for our manuscript beyond creating the `sce`
object we previously described.

The 10x Genomics file structure is replicated inside
[Histology](Analysis/Histology/) where we saved the image segmentation
analysis output to estimate the number of cells per spot. This involved
running some [histology image segmentation
software](https://www.mathworks.com/help/images/color-based-segmentation-using-k-means-clustering.html)
and the [counting code](Analysis/Histology/code) that requires our file
structure (`sgeID` input).

The main layer-level analysis code is located at
[Layer\_Guesses](Analysis/Layer_Guesses), for example
[layer\_specificity.R](Analysis/Layer_Guesses/layer_specificity.R) is
the R script for pseudo-bulking the spot-level data to create the
layer-level data. The `spatialLIBD` layer annotation files are saved in
the [First\_Round](Analysis/Layer_Guesses/First_Round/) and
[Second\_Round](Analysis/Layer_Guesses/Second_Round) directories which
you can upload to the shiny web application.

We also include directories with code for processing external datasets
such as [he\_layers](Analysis/he_layers),
[allen\_data](Analysis/allen_data/),
[hafner\_vglut](Analysis/hafner_vglut/).

We would like to highlight that a lot of the plotting code and
functionality from these scripts has been implemented in
*[spatialLIBD](https://bioconductor.org/packages/3.11/spatialLIBD)*
which would make a lot of our analysis simpler. Finally, for
reproducibility purposes we included the R session information in many
of our R scripts. Although in general we used R 3.6.1 and 3.6.2 with
Bioconductor release 3.10.

## [outputs](outputs/) directory

Contains outputs from the different unsupervised, semi-supervised, known
gene marker based and other clustering results. The analysis code that
generates these CSV files is located inside R Markdown files at the
[Analysis](Analysis/) directory such as
[SpatialDE\_clustering.Rmd](Analysis/SpatialDE_clustering.Rmd).

<a href="https://www.libd.org/"><img src="http://lcolladotor.github.io/img/LIBD_logo.jpg" width="250px"></a>

# Bibliography

\[1\] L. Collado-Torres, K. R. Maynard, and A. E. Jaffe. *LIBD Visium
spatial transcriptomics human pilot data inspector*.
<https://github.com/LieberInstitute/spatialLIBD> - R package version
1.0.0. 2020. DOI: 10.18129/B9.bioc.spatialLIBD. \<URL:
<http://www.bioconductor.org/packages/spatialLIBD>\>.

\[2\] K. R. Maynard, L. Collado-Torres, L. M. Weber, C. Uytingco, et al.
“Transcriptome-scale spatial gene expression in the human dorsolateral
prefrontal cortex”. In: *bioRxiv* (2020). DOI:
10.1101/2020.02.28.969931. \<URL:
<https://www.biorxiv.org/content/10.1101/2020.02.28.969931v1>\>.

# Internal

  - JHPCE location:
    `/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot`
  - Main `sce` R object file:
    `/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Human_DLPFC_Visium_processedData_sce_scran.Rdata`.

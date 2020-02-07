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

genes_to_plot <- c(
    'ATP1A2',
    'FABP7',
    'ADCYAP1',
    'NEFM',
    ## Originally it was NEFN, but that doesn't exist, fixed the typo
    'NEFH',
    'PVALB',
    'VAT1L',
    'NTNG2',
    'CUX2',
    'ENC1',
    'RORB',
    'RXFP1',
    'RPRM',
    'ETV1',
    'PCP4',
    'B3GALT2',
    'CCK',
    'MBP' ## From the figure, not on the Slack message
)

## Load sig_genes data
load('rda/layer_sig_genes.Rdata', verbose = TRUE)

## Looking for NEFN
sort(unique(sig_genes$gene[grep('ne', tolower(sig_genes$gene))]))
# [1] "NECAB1"   "NECAB2"   "NEFH"     "NEFL"     "NEFM"     "NEUROD6"  "SERPINE2"

## After editing genes_to_plot, we are good to go
stopifnot(all(genes_to_plot %in% sig_genes$gene))

## Prepare the data needed for making the plots
sig_genes_sub <- sig_genes[match(genes_to_plot, sig_genes$gene), ]
sig_genes_unique <- splitit(sig_genes_sub$ensembl)

## For the titles
sig_genes_df <- sig_genes_sub
sig_genes_df$in_rows <-
    sapply(sig_genes_df$in_rows, paste0, collapse = ';')
sig_genes_df$results <-
    sapply(sig_genes_df$results, paste0, collapse = ';')


## Customize plotting code from
# source(here('Analysis', 'spatialLIBD_global_plot_code.R'))
sce_image_clus_gene_p <-
    function(sce, d, sampleid, spatial, title) {
        p <-
            ggplot(d,
                aes(
                    x = imagecol,
                    y = imagerow,
                    fill = COUNT,
                    color = COUNT,
                    key =  key
                ))
        
        if (spatial) {
            p <-
                p + geom_spatial(
                    data = subset(metadata(sce)$image, sample == sampleid),
                    aes(grob = grob),
                    x = 0.5,
                    y = 0.5
                )
        }
        p <- p +
            geom_point(shape = 21,
                size = 1.25,
                stroke = 0.25) +
            # geom_point(size = 1.25, alpha = 0.25) +
            coord_cartesian(expand = FALSE) +
            scale_fill_gradientn(
                colors = viridis(21)
            ) +
            scale_color_gradientn(
                colors = viridis(21)
            ) +
            xlim(0, max(sce$width)) +
            ylim(max(sce$height), 0) +
            xlab("") + ylab("") +
            labs(fill = "COUNT") +
            ggtitle(title) +
            theme_set(theme_bw(base_size = 10)) +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black"),
                axis.text = element_blank(),
                axis.ticks = element_blank()
            )
        return(p)
    }



## Make gene grid plots
## Takes about 1 - 1.5 hours
pdf_dir <- 'pdf/gene_grid/RNAscope'
dir.create(pdf_dir, showWarnings = FALSE, recursive = TRUE)

## Only make the plots for the unique ones
## and only for the last 4 samples
samples_to_plot <- tail(unique(sce$sample_name), 4)
assayname <- 'logcounts'

for (j in samples_to_plot) {
    # j <-  samples_to_plot[1]
    dir.create(file.path(pdf_dir, j), showWarnings = FALSE)
    
    max_UMI <-
        max(assays(sce)[[assayname]][names(sig_genes_unique), sce$sample_name %in% j])
    max_UMI
    x <-
        assays(sce)[[assayname]][names(sig_genes_unique), sce$sample_name %in% j]
    min_UMI <- min(as.vector(x)[as.vector(x) > 0])
    min_UMI
    
    for (i in match(names(sig_genes_unique), sig_genes_sub$ensembl)) {
        # i <- 1
        # i <- 15 ## PCP4
        # i <- 11 ## RORB
        message(paste(
            Sys.time(),
            'making the plot for',
            i,
            'gene',
            sig_genes_sub$gene[i]
        ))
        
        p <- sce_image_grid_gene(
            sce[, sce$sample_name == j],
            geneid = paste0(sig_genes_sub$gene[i], '; ', sig_genes_sub$ensembl[i]),
            return_plots = TRUE,
            ... = gsub('top', 'r', gsub(
                'Layer', 'L', sig_genes_df$results[i]
            )),
            spatial = TRUE,
            assayname = assayname,
            minCount = 0
        )
        
        p2 <- p[[1]] + scale_fill_gradientn(
                colors = c('aquamarine4', 'springgreen', 'goldenrod', 'red'), na.value = add.alpha('black', 0.175), name = assayname, values = scales::rescale(c(min_UMI, 2, 4, max_UMI))
            ) +
            scale_color_gradientn(
                colors =  c('aquamarine4', 'springgreen', 'goldenrod', 'red'), na.value =  add.alpha('black', 0.175), name = assayname, values = scales::rescale(c(min_UMI, 2, 4, max_UMI))
            )
        pdf(
            file.path(
                pdf_dir,
                j,
                paste0(
                    sig_genes_sub$gene[i],
                    gsub('top', 'r', gsub(
                        'Layer', 'L', sig_genes_df$results[i]
                    )),
                    '.pdf'
                )
            ),
            useDingbats = FALSE,
            height = 8,
            width = 9
        )
        print(p2)
        dev.off()
    }
}

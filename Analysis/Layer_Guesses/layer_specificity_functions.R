## From https://gist.githubusercontent.com/mages/5339689/raw/2aaa482dfbbecbfcb726525a3d81661f9d802a8e/add.alpha.R
add.alpha <- function(col, alpha = 1) {
    if (missing(col))
        stop("Please provide a vector of colours.")
    apply(sapply(col, col2rgb) / 255, 2,
        function(x)
            rgb(x[1], x[2], x[3], alpha = alpha))
}

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
genes_RNAscope_raw <-
    read_xlsx(here('Analysis', 'RNAscope_Probe_List_December2018.xlsx'),
        sheet = 'Human Probes')

## Build gene annotation data.frame for heatmap
gene_ann <- function(x) {
    m_km <- match(tolower(x), tolower(genes_km_raw$Gene))
    m_bm <- match(tolower(x), tolower(genes_bm_raw$Gene))
    m_RNAscope <-
        match(tolower(x), tolower(genes_RNAscope_raw[['Gene Symbol']]))
    res <-
        data.frame(
            KM_Zeng = factor(!is.na(m_km), levels = c('FALSE', 'TRUE')),
            BM = factor(!is.na(m_bm), levels = c('FALSE', 'TRUE')),
            RNAscope = factor(!is.na(m_RNAscope), levels = c('FALSE', 'TRUE'))
        )
    rownames(res) <- make.unique(x)
    return(res)
}

ann_colors <-
    list(
        BM = c(
            'FALSE' = RColorBrewer::brewer.pal(6, 'Dark2')[1],
            'TRUE' = RColorBrewer::brewer.pal(6, 'Dark2')[2]
        ),
        KM_Zeng = c(
            'FALSE' = RColorBrewer::brewer.pal(6, 'Dark2')[3],
            'TRUE' = RColorBrewer::brewer.pal(6, 'Dark2')[4]
        ),
        RNAscope = c(
            'FALSE' = RColorBrewer::brewer.pal(6, 'Dark2')[5],
            'TRUE' = RColorBrewer::brewer.pal(6, 'Dark2')[6]
        )
    )

## Plotting code
plot_markers_logfc <-
    function(x,
        pval.type = c('any', 'all'),
        prefix = 'logFC',
        ...) {
        lapply(seq_along(x), function(chosen) {
            interesting <- x[[chosen]]
            if (pval.type == 'any') {
                best.set <-
                    interesting[interesting$Top <= 6, ] ## for pval.type = 'any'
            } else {
                best.set <- head(interesting, 30) ## for pval.type == 'all'
            }
            logFCs <- getMarkerEffects(best.set, prefix = prefix)
            print(
                pheatmap(
                    logFCs,
                    main = names(x)[chosen],
                    # color = colorRampPalette(c("white", "blue"))(100),
                    annotation_colors = ann_colors,
                    annotation_row = gene_ann(rownames(logFCs)),
                    annotation_names_row = TRUE,
                    ...
                )
            )
            return(NULL)
        })
    }



## Copy the sce object but with other
## rownames so it'll be easier to re-use existing functions
sce_layer_symbol <- sce_layer
rownames(sce_layer_symbol) <-
    make.unique(rowData(sce_layer)$gene_name)

plot_markers_expr <-
    function(x,
        pval.type = c('any', 'all'),
        prefix = NULL,
        ...) {
        lapply(seq_along(x), function(chosen) {
            interesting <- x[[chosen]]
            if (pval.type == 'any') {
                best.set <-
                    interesting[interesting$Top <= 6,] ## for pval.type = 'any'
            } else {
                best.set <- head(interesting, 30) ## for pval.type == 'all'
            }
            pheat <- plotHeatmap(
                sce_layer_symbol,
                features = rownames(best.set),
                main = names(x)[chosen],
                colour_columns_by = c(
                    'subject',
                    'subject_position',
                    'sample_name',
                    'c_k5_k7',
                    'c_k7_k7',
                    'c_k20_k7',
                    'kmeans_k7',
                    'layer_guess'
                ),
                # color = colorRampPalette(c("white", "blue"))(100),
                # annotation_colors = ann_colors,
                # color = viridis::viridis(21),
                annotation_row = gene_ann(rownames(best.set)),
                annotation_names_row = TRUE,
                ...
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

plot_markers_loop <-
    function(x,
        pdf_header,
        FUN,
        h = 14,
        prefix = NULL,
        ...) {
        for (pval in names(x)) {
            for (direc in names(x[[1]])) {
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
                    height = h
                )
                FUN(x[[pval]][[direc]], pval.type = pval, prefix = prefix, ...)
                dev.off()
            }
        }
    }

find_marker_gene <-
    function(x,
        markers,
        pval = 'any',
        direc = 'any',
        layer = 1) {
        m <- match(x, rownames(markers[[pval]][[direc]][[layer]]))
        res <- markers[[pval]][[direc]][[layer]][m,]
        res$rownum <- m
        return(res)
    }


## Write a function for extracting the data
sig_genes_extract <- function(tstats, pvals, n = 10) {
    stopifnot(identical(dim(tstats), dim(pvals)))
    sig_genes <- apply(tstats, 2, function(x) {
        rowData(sce_layer)$gene_name[order(x, decreasing = TRUE)[1:10]]
    })

    sig_i <- apply(tstats, 2, function(x) {
        order(x, decreasing = TRUE)[seq_len(n)]
    })
    sig_genes_tstats <-
        sapply(seq_len(ncol(sig_i)), function(i) {
            tstats[sig_i[, i], i]
        })
    sig_genes_pvals <-
        sapply(seq_len(ncol(sig_i)), function(i) {
            pvals[sig_i[, i], i]
        })
    sig_genes_fdr <-
        sapply(seq_len(ncol(sig_i)), function(i) {
            apply(pvals, 2, p.adjust, 'fdr')[sig_i[, i], i]
        })
    dimnames(sig_genes_fdr) <-
        dimnames(sig_genes_tstats) <-
        dimnames(sig_genes_pvals) <- dimnames(sig_genes)

    ## Combine into a long format table
    sig_genes_tab <- data.frame(
        top = rep(seq_len(n), n = ncol(tstats)),
        layer = rep(colnames(sig_genes), each = n),
        gene = as.character(sig_genes),
        tstat = as.numeric(sig_genes_tstats),
        pval = as.numeric(sig_genes_pvals),
        fdr = as.numeric(sig_genes_fdr),
        gene_index = as.integer(sig_i),
        stringsAsFactors = FALSE
    )
    sig_genes_tab$ensembl <-
        rownames(sce_layer)[sig_genes_tab$gene_index]

    ## Add gene marker labels
    sig_genes_tab <-
        cbind(sig_genes_tab, gene_ann(sig_genes_tab$gene))
    rownames(sig_genes_tab) <- NULL
    return(sig_genes_tab)
}

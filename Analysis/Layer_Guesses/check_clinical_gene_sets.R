###
library('readxl')
library('limma')
library('sessioninfo')
library('parallel')
library('jaffelab')
library('janitor')
library('lattice')
library('org.Hs.eg.db')
library('GenomicFeatures')
library('scran')
library('here')

## load sce object
sce_layer_file <-
    here('Analysis', 'Layer_Guesses', 'rda', 'sce_layer.Rdata')
if (file.exists(sce_layer_file))
    load(sce_layer_file, verbose = TRUE)
	
###################
## load modeling outputs
load("rda/eb_contrasts.Rdata")
load("rda/eb0_list.Rdata")

## Extract the p-values
logFC0_contrasts <- sapply(eb0_list, function(x) {
    x$coef[, 2, drop = FALSE]
})
rownames(logFC0_contrasts) = rownames(eb_contrasts)
## Extract the p-values
pvals0_contrasts <- sapply(eb0_list, function(x) {
    x$p.value[, 2, drop = FALSE]
})
rownames(pvals0_contrasts) = rownames(eb_contrasts)
fdrs0_contrasts = apply(pvals0_contrasts, 2, p.adjust, "fdr")

## Extract the t-stats
t0_contrasts <- sapply(eb0_list, function(x) {
    x$t[, 2, drop = FALSE]
})
rownames(t0_contrasts) = rownames(eb_contrasts)

#########################
## load in gene sets ####
#########################

##################################
## Satterstrom et al, Cell 2020 ##
##################################
asd_exome = read_excel("gene_sets/1-s2.0-S0092867419313984-mmc2.xlsx", 
	sheet = 2)
asd_exome = as.data.frame(asd_exome)

## get ensembl IDs
asd_exome_geneList = apply(asd_exome[,
    c("ASC33_2014",
        "SSC27_2014",
        "ASC65_2015",
        "ASC102_2018",
        "ASD53",
        "DDID49")], 2,
    function(x)
        asd_exome$ensembl_gene_id[x == 1])
names(asd_exome_geneList) = gsub("_", ".", names(asd_exome_geneList))
names(asd_exome_geneList) = paste0("Gene_Satterstrom_",
    names(asd_exome_geneList))

###############
### SFARI #####
###############

asd_sfari = read.csv("gene_sets/SFARI-Gene_genes_01-03-2020release_02-04-2020export.csv",
    as.is = TRUE)
asd_sfari_geneList = list(
    Gene_SFARI_all = asd_sfari$ensembl.id,
    Gene_SFARI_high = asd_sfari$ensembl.id[asd_sfari$gene.score < 3],
    Gene_SFARI_syndromic = asd_sfari$ensembl.id[asd_sfari$syndromic == 1]
)

# #################
# ## harmonizome ##
# #################

# harmonizome = read.delim(
    # "gene_sets/Harmonizome_CTD Gene-Disease Associations Dataset.txt",
    # as.is = TRUE,
    # skip = 1
# )
# ## add ensembl
# ens = select(org.Hs.eg.db,
    # columns = c("ENSEMBL", "ENTREZID"),
    # keys = as.character(unique(harmonizome$GeneID)))
# harmonizome$ensemblID = ens$ENSEMBL[match(harmonizome$GeneID, ens$ENTREZID)]

# ## split by dx
# harmonizome_geneList = split(harmonizome$ensemblID, harmonizome$Disease)

# ## filter by set size
# harmonizome_geneList = harmonizome_geneList[lengths(harmonizome_geneList) >= 100]
# names(harmonizome_geneList) = gsub(" ", ".", names(harmonizome_geneList))
# names(harmonizome_geneList) = paste0("Harmonizome_",
    # names(harmonizome_geneList))

####################
### birnbaum sets ##
####################

birnbaum = read_excel("gene_sets/Supplementary Tables for paper.Birnbaum November 2013.AJP.xlsx",
    sheet = 1)
ens2 = select(org.Hs.eg.db,
    columns = c("ENSEMBL", "ENTREZID"),
    keys = as.character(unique(birnbaum$`EntrezGene ID`)))
birnbaum$ensemblID = ens2$ENSEMBL[match(birnbaum$`EntrezGene ID`, ens2$ENTREZID)]

birnbaum_geneList = split(birnbaum$ensemblID, birnbaum$`Gene Set`)
names(birnbaum_geneList) = gsub(" ", ".", names(birnbaum_geneList))
names(birnbaum_geneList) = gsub("-", ".", names(birnbaum_geneList))
names(birnbaum_geneList) = paste0("Gene_Birnbaum_",
    names(birnbaum_geneList))
birnbaum_geneList = birnbaum_geneList[rev(seq(along=birnbaum_geneList))]

######################
## psychENCODE DEGs ##
######################

psychENCODE = as.data.frame(read_excel("gene_sets/aat8127_Table_S1.xlsx", sheet = "DGE"))

pe_geneList = with(
    psychENCODE,
    list(
        DE_PE_ASD.Up = ensembl_gene_id[ASD.t.value > 0 & ASD.fdr < 0.05],
        DE_PE_ASD.Down = ensembl_gene_id[ASD.t.value < 0 & ASD.fdr < 0.05],
        DE_PE_BD.Up = ensembl_gene_id[BD.t.value > 0 & BD.fdr < 0.05],
        DE_PE_BD.Down = ensembl_gene_id[BD.t.value < 0 & BD.fdr < 0.05],
		DE_PE_SCZ.Up = ensembl_gene_id[SCZ.t.value > 0 & SCZ.fdr < 0.05],
        DE_PE_SCZ.Down = ensembl_gene_id[SCZ.t.value < 0 & SCZ.fdr < 0.05]
		)
)

#################
## brainseq  ####
#################

## DLPFC RiboZero
load(
    "/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_noHGoldQSV_matchDLPFC.rda"
)

bs2_geneList = with(outGene,
    list(DE_BS2_SCZ.Up = ensemblID[logFC > 0 & adj.P.Val < 0.05],
        DE_BS2_SCZ.Down = ensemblID[logFC < 0 & adj.P.Val < 0.05]))


##############################
### Sestan DS Neuron 2017? ###

ds = read_excel("gene_sets/1-s2.0-S0896627316000891-mmc4.xlsx",skip=2)
ds = clean_names(ds)
ds = as.data.frame(ds)
ens3 = select(org.Hs.eg.db,
    columns = c("ENSEMBL", "ENTREZID"),
    keys = as.character(unique(ds$geneid)))
ds$ensemblID = ens3$ENSEMBL[match(ds$geneid, ens3$ENTREZID)]
ds$fold_difference_log2 = as.numeric(ds$fold_difference_log2)
ds$p_value = readr::parse_number(ds$p_value)
ds$qval = readr::parse_number(ds$qval)

ds_geneList = list(DE_DS_DS.Up = ds$ensemblID[ds$fold_difference_log2 > 0 & ds$qval < 0.05],
        DE_DS_DS.Down = ds$ensemblID[ds$fold_difference_log2 < 0 & ds$qval < 0.05])

#############################
## various TWAS sets ########
#############################

## brainseq 2
load("/dcl01/ajaffe/data/lab/dg_hippo_paper/rdas/tt_objects_gene.Rdata")
tt_dlpfc=  as.data.frame(tt[tt$region == "DLPFC",])
tt_dlpfc$ensemblID = ss(tt_dlpfc$geneid, "\\.")

## PE 
twas_sczd = as.data.frame(read_excel("gene_sets/aat8127_Table_S4.xlsx", sheet = "SCZ.TWAS"))
twas_sczd$TWAS.FDR = p.adjust(twas_sczd$TWAS.P, "fdr")
twas_asd = as.data.frame(read_excel("gene_sets/aat8127_Table_S4.xlsx", sheet = "ASD.TWAS"))
twas_asd$TWAS.FDR = p.adjust(twas_asd$TWAS.P, "fdr")
twas_bpdscz = as.data.frame(read_excel("gene_sets/aat8127_Table_S4.xlsx", sheet = "BD.SCZ"))
twas_bpdscz$TWAS.FDR = p.adjust(twas_bpdscz$TWAS.P, "fdr")

twas_geneList = list(TWAS_BS2_SCZ.Up = tt_dlpfc$ensemblID[tt_dlpfc$TWAS.Z > 0 & tt_dlpfc$TWAS.FDR < 0.05],
			TWAS_BS2_SCZ.Down = tt_dlpfc$ensemblID[tt_dlpfc$TWAS.Z < 0 & tt_dlpfc$TWAS.FDR < 0.05],
			TWAS_PE_SCZ.Up = twas_sczd$GeneID[twas_sczd$TWAS.Z > 0 & twas_sczd$TWAS.FDR < 0.05],
			TWAS_PE_SCZ.Down = twas_sczd$GeneID[twas_sczd$TWAS.Z < 0 & twas_sczd$TWAS.FDR < 0.05],
			TWAS_PE_ASD.Up = twas_asd$GeneID[twas_asd$TWAS.Z > 0 & twas_asd$TWAS.FDR < 0.05],
			TWAS_PE_ASD.Down = twas_asd$GeneID[twas_asd$TWAS.Z < 0 & twas_asd$TWAS.FDR < 0.05],
			TWAS_PE_SCZBD.Up = twas_bpdscz$ID[twas_bpdscz$TWAS.Z > 0 & twas_bpdscz$TWAS.FDR < 0.05],
			TWAS_PE_SCZBD.Down = twas_bpdscz$ID[twas_bpdscz$TWAS.Z < 0 & twas_bpdscz$TWAS.FDR < 0.05])



###############
### combine ###
###############

## gene list ##
geneList = c(
    birnbaum_geneList,
	asd_sfari_geneList,
	asd_exome_geneList,
    pe_geneList,
    bs2_geneList,
    ds_geneList,
	twas_geneList
)

## filter for those present in spatial data
geneList_present = lapply(geneList, function(x) {
    x = x[!is.na(x)]
    x[x %in% rownames(t0_contrasts)]
})

## do enrichment
enrich_stat_list = eb0_list
for (i in seq(along = eb0_list)) {
    layer = t0_contrasts[, i] > 0 & fdrs0_contrasts[, i] < 0.1
	tabList = mclapply(geneList_present, function(g) {
        tt = table(Set = factor(names(layer) %in% g, c(FALSE, TRUE)),
            Layer = factor(layer, c(FALSE, TRUE)))
    }, mc.cores = 8)
	enrichList = lapply(tabList,fisher.test)
	
    o = data.frame(
        OR = sapply(enrichList, "[[", "estimate"),
        Pval = sapply(enrichList, "[[", "p.value"),
		NumSig = sapply(tabList, function(x) x[2,2])
    )
    rownames(o) = gsub(".odds ratio", "", rownames(o))
    enrich_stat_list[[i]] = o
}
enrichTab = do.call("cbind", enrich_stat_list)

#  name
enrichTab$Type = ss(rownames(enrichTab), "_", 1)
enrichTab$Type[enrichTab$Group == "Birnbaum"] = "Birnbaum"
enrichTab$Type[enrichTab$Type == "Gene"] = "ASD"
enrichTab$Group = ss(rownames(enrichTab), "_", 2)
enrichTab$Set = ss(rownames(enrichTab), "_", 3)
enrichTab$ID = rownames(enrichTab)
enrichTab$SetSize = sapply(geneList_present, length)

## look at enrichment
pMat = enrichTab[, grep("Pval", colnames(enrichTab))]
orMat = enrichTab[, grep("OR", colnames(enrichTab))]
colnames(pMat) = ss(colnames(pMat), "\\.")
colnames(orMat) = ss(colnames(orMat), "\\.")
pMat < 0.05 / nrow(pMat)
pMat < 0.001
round(-log10(pMat),1)

######################
## pull out results ##
######################

## summary stats from genes
enrichTab["Gene_SFARI_all",]
enrichTab["Gene_Satterstrom_ASC102.2018",]
enrichTab["Gene_Satterstrom_ASD53",]
enrichTab["Gene_Satterstrom_DDID49",]

## Satterstrom deep dive
sat_102_l2= which(t0_contrasts[,"Layer2"] > 0 & fdrs0_contrasts[,"Layer2"] < 0.1 & 
		rownames(t0_contrasts) %in% geneList_present$Gene_Satterstrom_ASC102.2018)
rowData(sce_layer)$gene_name[sat_102_l2]
sat_102_l5= which(t0_contrasts[,"Layer5"] > 0 & fdrs0_contrasts[,"Layer5"] < 0.1 & 
		rownames(t0_contrasts) %in% geneList_present$Gene_Satterstrom_ASC102.2018)
rowData(sce_layer)$gene_name[sat_102_l5]

sat_49_l2= which(t0_contrasts[,"Layer2"] > 0 & fdrs0_contrasts[,"Layer2"] < 0.1 & 
		rownames(t0_contrasts) %in% geneList_present$Gene_Satterstrom_DDID49)
cat(rowData(sce_layer)$gene_name[sat_49_l2], sep=", ")

sat_53_l5= which(t0_contrasts[,"Layer5"] > 0 & fdrs0_contrasts[,"Layer5"] < 0.1 & 
		rownames(t0_contrasts) %in% geneList_present$Gene_Satterstrom_ASD53)
cat(rowData(sce_layer)$gene_name[sat_53_l5], sep=", ")

## case control - asd
enrichTab["DE_PE_ASD.Up",]
enrichTab["DE_PE_ASD.Down",]

## case control - sczd
enrichTab[c("DE_PE_SCZ.Up","DE_BS2_SCZ.Up"),]
enrichTab[c("DE_PE_SCZ.Down","DE_BS2_SCZ.Down"),]


#############################
### GSEA #####################
#############################

## do enrichment
gst_tab = apply(t0_contrasts, 2, function(tt) {
	sapply(geneList_present, function(g) {
       geneSetTest(index = rownames(t0_contrasts) %in% g,
		statistics = tt, alternative = "up")
    })
})
round(-log10(gst_tab),1)

## check densities 
mypar(ncol(t0_contrasts),1)
g = rownames(t0_contrasts) %in% geneList_present$Gene_Satterstrom_ASC102.2018 
g_asd = rownames(t0_contrasts) %in% geneList_present$Gene_Satterstrom_ASD53
g_dd = rownames(t0_contrasts) %in% geneList_present$Gene_Satterstrom_DDID49
for(i in 1:ncol(t0_contrasts)) {
	plot(density(t0_contrasts[!g,i]), lwd=2,col="black",xlab="",
		main=colnames(t0_contrasts)[i],xlim = c(-10,15))
	lines(density(t0_contrasts[g,i]), lwd=2,col="red")
	lines(density(t0_contrasts[g_asd,i]), lwd=2,col="red",lty=2)
	lines(density(t0_contrasts[g_dd,i]), lwd=2,col="red",lty=3)
	abline(v=0,lty=2)
}

pdf("pdf/ASD_genes_layer_density.pdf",h=4,useDingbats=FALSE)
par(mar=c(5,6,1,1),cex.axis= 1.4,cex.lab=1.8)
for(i in 1:ncol(t0_contrasts)) {
    layer = t0_contrasts[, i] > 0 & fdrs0_contrasts[, i] < 0.1
	plot(density(t0_contrasts[!g,i]), lwd=3,col="black",
		xlab=paste0(colnames(t0_contrasts)[i], ": Specificity T-stats"),
		sub = "", main="",xlim = c(-8,8))
	lines(density(t0_contrasts[g,i]), lwd=3,col="red")
	lines(density(t0_contrasts[g_asd,i]), lwd=3,col="red",lty=2)
	lines(density(t0_contrasts[g_dd,i]), lwd=3,col="red",lty=3)
	abline(v=0,lty=2)
	abline(v=	min(t0_contrasts[layer,i]))

	ll = ifelse(i == 1, "topright", "topleft")
	legend(ll, c("BG", "102 All", "53 ASD", "49 DDID"), bty="n",
		col = c("black","red","red","red"),	lty = c(1,1,2,3),cex=1.5,lwd=4)
}
dev.off()

diag(cor(t(-log10(gst_tab)),t(-log10(pMat))))
################
## make plots ##
################

################
## dotplots ####
################

library(ggplot2)

## make long
enrichLong = reshape2::melt(enrichTab[,c(1,3,5,7,9,11,13,15:18)],id.vars = 8:11)
colnames(enrichLong)[5:6] = c("Layer", "OR")
enrichLong_P = reshape2::melt(enrichTab[,c(2,4,6,8,10,12,14,15:18)],id.vars = 8:11)
identical(enrichLong$ID, enrichLong_P$ID)
enrichLong$P = enrichLong_P$value
enrichLong$Layer = ss(as.character(enrichLong$Layer), "\\.")
enrichLong$ID = factor(enrichLong$ID, levels=rev(rownames(enrichTab)))
enrichLong$Set = factor(enrichLong$Set, levels=unique(rev(enrichTab$Set)))

## overall ##
enrichLong$P_thresh = enrichLong$P
enrichLong$P_thresh[enrichLong$P_thresh < 2.2e-16] = 2.2e-16

#### overall ###
pdf("pdf/clinical_gene_sets_all_dotplot.pdf",h=9.5,w=7.5)
print(
	ggplot(enrichLong, aes(x=Layer, y=ID, size=OR, color=-log10(P_thresh))) +
		geom_point() + scale_color_continuous(low="white", high="darkred", 
			name = "-log10(p)", guide=guide_colorbar(reverse=TRUE))+ 
			ylab(NULL) +xlab(NULL) + ggtitle("Clinical Gene Enrichment") + 
			theme_dark() + 
			theme(text = element_text(size = 20),
				axis.text.x = element_text(angle = 90, hjust = 1),
				legend.key = element_rect(colour = "transparent", fill = "white")) +	
			# guides(size = guide_legend(override.aes = list(color = "white"))) + 				
			scale_size(range=c(1, 10))  
)
dev.off()		

#### ASD focus ###
enrichLong_typeList = split(enrichLong, enrichLong$Type)
enrichLong_typeList = lapply(enrichLong_typeList, function(x) {
	x$ID = paste0(x$Group, ":", x$Set)
	x$ID = factor(x$ID, unique(rev(paste0(enrichTab$Group, ":", enrichTab$Set))))
	return(x)
})

pdf("pdf/clinical_gene_sets_byType_dotplot.pdf",h=6,w=7.5)
for(i in seq(along=enrichLong_typeList)) {
	print(
		ggplot(enrichLong_typeList[[i]], aes(x=Layer, y=ID, size=OR, color=-log10(P_thresh))) +
			geom_point() + scale_color_continuous(low="white", high="darkred", 
			name = "-log10(p)", guide=guide_colorbar(reverse=TRUE))+ 
			ylab(NULL) +xlab(NULL) + ggtitle(paste(names(enrichLong_typeList)[i],
				"Enrichment")) + 
			theme_dark() + 
			theme(text = element_text(size = 20),
				axis.text.x = element_text(angle = 90, hjust = 1),
				legend.key = element_rect(colour = "transparent", fill = "white")) +	
			scale_size(range=c(1, 10))
		)
}		
dev.off()

enrichLong_groupList = split(enrichLong, paste(enrichLong$Type, enrichLong$Group))
names(enrichLong_groupList)[3] = "Birnbaum Gene"

pdf("pdf/clinical_gene_sets_byGroup_dotplot.pdf",h=6,w=7.5)
for(i in seq(along=enrichLong_groupList)) {
	print(
		ggplot(enrichLong_groupList[[i]],
			aes(x=Layer, y=Set, size=OR, color=-log10(P_thresh))) +
			geom_point() + scale_color_continuous(low="white", high="darkred", 
			name = "-log10(p)", guide=guide_colorbar(reverse=TRUE))+ 
			ylab(NULL) +xlab(NULL) + ggtitle(paste(names(enrichLong_groupList)[i],
				"Enrichment")) + 
			theme_dark() + 
			theme(text = element_text(size = 20),
				axis.text.x = element_text(angle = 90, hjust = 1),
				legend.key = element_rect(colour = "transparent", fill = "white")) +	
			scale_size(range=c(1, 10)) + guides(fill = "white") 
		)
}		
dev.off()



## ggplot2 style?
negp_long = reshape2::melt(-log10(pMat))

theSeq = seq(0, 12, by = 0.1)
my.col <- colorRampPalette(c("white", "red"))(length(theSeq))

## overall

negPmat[negPmat > 12] = 12

## ASD
negPmat_ASD = negPmat[c(7:9, 4:6, 24:25), ]
rownames(negPmat_ASD)[7:8] = c("DE_ASD_Up", "DE_ASD_Down")

negPmat_ASD = as.matrix(negPmat_ASD)
negPmat_ASD = negPmat_ASD[, c(2:7, 1)]
negPmat_ASD = negPmat_ASD[nrow(negPmat_ASD):1, ]

pdf("pdf/ASD_risk_gene_heatmap.pdf", width = 8)
print(
    levelplot(
        t(negPmat_ASD),
        asPEct = "fill",
        at = theSeq,
        col.regions = my.col,
        ylab = "",
        xlab = "",
        scales = list(x = list(rot = 90, cex = 1.5), y = list(cex = 1.5))
    )
)
dev.off()

enrichTab_harm = enrichTab[enrichTab$Group == "Harmonizome", ]
fdrTab_harm = apply(enrichTab_harm[, grep("Pval", colnames(enrichTab_harm))], 2, p.adjust, "fdr")
colSums(fdrTab_harm < 0.05)
colSums(enrichTab_harm[, grep("Pval", colnames(enrichTab_harm))] < 1e-10)
colSums(enrichTab_harm[, grep("Pval", colnames(enrichTab_harm))] < 1e-20)
harm_top10 = apply(enrichTab_harm[, grep("Pval", colnames(enrichTab_harm))], 2,
    function(x)
        gsub("Harmonizome_", "", rownames(enrichTab_harm)[order(x)[1:10]]))

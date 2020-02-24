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
library('RColorBrewer')
library('ggplot2')
library('fields')

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
enrichLong = reshape2::melt(enrichTab[,c(seq(1,19,by=3),22:26)],id.vars = 8:12)
colnames(enrichLong)[6:7] = c("Layer", "OR")
enrichLong_P = reshape2::melt(enrichTab[,c(seq(2,20,by=3),22:26)],id.vars = 8:12)
identical(enrichLong$ID, enrichLong_P$ID)
enrichLong$P = enrichLong_P$value
enrichLong$Layer = ss(as.character(enrichLong$Layer), "\\.")
enrichLong$ID = factor(enrichLong$ID, levels=rev(rownames(enrichTab)))
enrichLong$Set = factor(enrichLong$Set, levels=unique(rev(enrichTab$Set)))
enrichLong$FDR = p.adjust(enrichLong$P, "fdr")

## what p-value controls FDR?
enrichLongSort = enrichLong[order(enrichLong$P),]
max(enrichLongSort$P[enrichLongSort$FDR < 0.05] )
# 0.01009034

## overall ##
enrichLong$P_thresh = enrichLong$P
enrichLong$P_thresh[enrichLong$P_thresh < 2.2e-16] = 2.2e-16

### ASD focus
enrichLong_ASD = enrichLong[enrichLong$ID %in% 
	c("Gene_SFARI_all", "Gene_Satterstrom_ASC102.2018",
	"Gene_Satterstrom_ASD53", "Gene_Satterstrom_DDID49",
	"DE_PE_ASD.Down", "DE_PE_ASD.Up",
	"TWAS_PE_ASD.Up", "TWAS_PE_ASD.Down"),]
enrichLong_ASD$ID2 =  as.character(droplevels(enrichLong_ASD$Set))
enrichLong_ASD$ID2[enrichLong_ASD$ID2 == "all"] = "SFARI"
enrichLong_ASD$ID2[enrichLong_ASD$ID2 == "ASC102.2018"] = "ASC102"
enrichLong_ASD$ID2[enrichLong_ASD$ID == "DE_PE_ASD.Up"] = "DE.Up"
enrichLong_ASD$ID2[enrichLong_ASD$ID == "DE_PE_ASD.Down"] = "DE.Down"
enrichLong_ASD$ID2[enrichLong_ASD$ID == "TWAS_PE_ASD.Up"] = "TWAS.Up"
enrichLong_ASD$ID2[enrichLong_ASD$ID == "TWAS_PE_ASD.Down"] = "TWAS.Down"
enrichLong_ASD$ID2 = factor(enrichLong_ASD$ID2, unique(enrichLong_ASD$ID2))

enrichLong_ASD$LayerFac = factor(as.character(enrichLong_ASD$Layer), 
	c("WM", paste0("Layer", 6:1)))
enrichLong_ASD = enrichLong_ASD[order(enrichLong_ASD$ID2, enrichLong_ASD$LayerFac),]

#### overall ###
pdf("pdf/clinical_gene_sets_ASDfocus_dotplot.pdf",h=9.5,w=7.5)
print(
	ggplot(enrichLong_ASD, aes(y=LayerFac, x=ID2, size=OR, color=-log10(P_thresh))) +
		geom_point() + scale_color_continuous(low="white", high="darkred", 
			name = "-log10(p)", guide=guide_colorbar(reverse=TRUE))+ 
			ylab(NULL) +xlab(NULL) + ggtitle("ASD Gene Enrichment") + 
			theme_dark() + 
			theme(text = element_text(size = 20),
				axis.text.x = element_text(angle = 90, hjust = 1),
				legend.key = element_rect(colour = "transparent", fill = "white")) +	
			# guides(size = guide_legend(override.aes = list(color = "white"))) + 				
			scale_size(range=c(1, 10))  
)
dev.off()		


### custom heatmap

midpoint = function(x) x[-length(x)] + diff(x)/2

customLayerEnrichment = function(enrichTab , groups, xlabs, 
	Pthresh = 12, ORcut = 3, enrichOnly = FALSE,
	layerHeights = c(0,40,55,75,85,110,120,135),
	mypal = c("white", colorRampPalette(brewer.pal(9,"BuGn"))(50))) {

	wide_p = -log10( enrichTab[groups,grep("Pval", colnames(enrichTab))])
	wide_p[wide_p > Pthresh] = Pthresh
	wide_p = t(round(wide_p[,
		c("WM.Pval", "Layer6.Pval", "Layer5.Pval", "Layer4.Pval", "Layer3.Pval","Layer2.Pval", "Layer1.Pval")],2))

	wide_or = enrichTab[groups,grep("OR", colnames(enrichTab))]
	wide_or= round(t(wide_or[,
		c("WM.OR", "Layer6.OR", "Layer5.OR", "Layer4.OR", "Layer3.OR", "Layer2.OR", "Layer1.OR")]),1)
	if(enrichOnly) wide_p[wide_or < 1] = 0
	wide_or[wide_p < ORcut] = ""

	image.plot(x = seq(0,ncol(wide_p),by=1), y = layerHeights, z = as.matrix(t(wide_p)),
		col = mypal,xaxt="n", yaxt="n",xlab = "", ylab="")
	axis(2, c("WM", paste0("L", 6:1)), at = midpoint(layerHeights),las=1)
	axis(1, rep("", ncol(wide_p)), at = seq(0.5,ncol(wide_p)-0.5))
	text(x = seq(0.5,ncol(wide_p)-0.5), y=-1*max(nchar(xlabs))/2, xlabs,
		xpd=TRUE, srt=45,cex=2,adj= 1)
	abline(h=layerHeights,v=0:ncol(wide_p))
	text(x = rep(seq(0.5,ncol(wide_p)-0.5),each = nrow(wide_p)), 
		y = rep(midpoint(layerHeights), ncol(wide_p)),
		as.character(wide_or),cex=1.5,font=2)
}
	
pdf("pdf/asd_geneSet_heatmap.pdf",w=6)
par(mar=c(8,4.5,2.5,1), cex.axis=2,cex.lab=2)
groups = unique(as.character(enrichLong_ASD$ID))[1:6]
xlabs  = as.character(enrichLong_ASD$ID2[match(groups, enrichLong_ASD$ID)])
customLayerEnrichment(enrichTab, groups,xlabs, enrichOnly=TRUE)
abline(v=4,lwd=3)
text(x = 3, y = 142, c("ASD"), xpd=TRUE,cex=2.5,font=2)

dev.off()


pdf("pdf/sczd_geneSet_heatmap.pdf",w=8)
par(mar=c(8,4.5,2.5,1), cex.axis=2,cex.lab=2)

groups =c("DE_PE_SCZ.Up", "DE_PE_SCZ.Down", 
	"DE_BS2_SCZ.Up", "DE_BS2_SCZ.Down", 
	"TWAS_BS2_SCZ.Up", "TWAS_BS2_SCZ.Down", "TWAS_PE_SCZ.Up",
	"TWAS_PE_SCZ.Down")
xlabs = ss(gsub("_SCZ", "", groups), "_", 2)
customLayerEnrichment(enrichTab, groups,xlabs, enrichOnly=TRUE)
abline(v=4,lwd=3)
text(x = c(2,6), y = 142, c("SCZD-DE", "SCZD-TWAS"), xpd=TRUE,cex=2.5,font=2)
dev.off()














#######################
## old ggplot 2 code ##
#######################

enrichLong_ASD$P_thresh[enrichLong_ASD$P_thresh < 1e-12] = 1e-12
ggplot(enrichLong_ASD, aes(y=LayerFac, x=ID2, fill =-log10(P_thresh))) +
		geom_tile(height = rep(c(5,3,4,1,4,1,2),times=8)) + 
		scale_fill_gradient(low="white", high="darkred", 
			name = "-log10(p)", guide=guide_colorbar(reverse=TRUE))+ 
			ylab(NULL) +xlab(NULL) + ggtitle("ASD Gene Enrichment") + 
			theme(text = element_text(size = 20),
				axis.text.x = element_text(angle = 90, hjust = 1),
				legend.key = element_rect(colour = "transparent", fill = "white")) +	
			# guides(size = guide_legend(override.aes = list(color = "white"))) + 				
			scale_size(range=c(1, 10))  
dev.off()


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



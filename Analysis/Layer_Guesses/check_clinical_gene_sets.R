###
library('readxl')
library('limma')
library('sessioninfo')
library('parallel')
library('jaffelab')

library(lattice)

###################
## load modeling outputs
load("rda/eb_contrasts.Rdata")
load("rda/eb0_list.Rdata")

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
asd_exome = read_excel("gene_sets/1-s2.0-S0092867419313984-mmc2.xlsx",sheet=2)
asd_exome = as.data.frame(asd_exome)

## get ensembl IDs
asd_exome_geneList = apply(asd_exome[,
	c("ASC33_2014", "SSC27_2014", "ASC65_2015",
	"ASC102_2018", "ASD53", "DDID49")], 2, 
		function(x) asd_exome$ensembl_gene_id[x==1])
names(asd_exome_geneList) = gsub("_", ".", names(asd_exome_geneList))
names(asd_exome_geneList) = paste0("Satterstrom_", 
		names(asd_exome_geneList))		

###############
### SFARI #####
###############

asd_sfari = read.csv("gene_sets/SFARI-Gene_genes_01-03-2020release_02-04-2020export.csv",
				as.is=TRUE)
asd_sfari_geneList = list(SFARI_all = asd_sfari$ensembl.id,
		SFARI_high = asd_sfari$ensembl.id[asd_sfari$gene.score < 3],
		SFARI_syndromic = asd_sfari$ensembl.id[asd_sfari$syndromic == 1])
		
#################
## harmonizome ##
#################
library(org.Hs.eg.db)
library(GenomicFeatures)

harmonizome = read.delim("gene_sets/Harmonizome_CTD Gene-Disease Associations Dataset.txt",
	as.is=TRUE,skip=1)
## add ensembl
ens = select(org.Hs.eg.db, columns = c("ENSEMBL", "ENTREZID"),
	keys = as.character(unique(harmonizome$GeneID)))
harmonizome$ensemblID = ens$ENSEMBL[match(harmonizome$GeneID, ens$ENTREZID)]

## split by dx
harmonizome_geneList = split(harmonizome$ensemblID, harmonizome$Disease)

## filter by set size
harmonizome_geneList = harmonizome_geneList[lengths(harmonizome_geneList) >= 100]
names(harmonizome_geneList) = gsub(" ", ".", names(harmonizome_geneList))
names(harmonizome_geneList) = paste0("Harmonizome_",
	names(harmonizome_geneList)) 

####################
### birnbaum sets ##
####################

birnbaum = read_excel("gene_sets/Supplementary Tables for paper.Birnbaum November 2013.AJP.xlsx",sheet=1)
ens2 = select(org.Hs.eg.db, columns = c("ENSEMBL", "ENTREZID"),
	keys = as.character(unique(birnbaum$`EntrezGene ID`)))
birnbaum$ensemblID = ens2$ENSEMBL[match(birnbaum$`EntrezGene ID`, ens2$ENTREZID)]

birnbaum_geneList = split(birnbaum$ensemblID, birnbaum$`Gene Set`)
names(birnbaum_geneList) = gsub(" ", ".", names(birnbaum_geneList))
names(birnbaum_geneList) = gsub("-", ".", names(birnbaum_geneList))
names(birnbaum_geneList) = paste0("Birnbaum_",
	names(birnbaum_geneList)) 

######################
## psychENCODE DEGs ##
######################

psychENCODE = as.data.frame(read_excel("gene_sets/aat8127_Table_S1.xlsx", sheet = "DGE"))
rownames(psychENCODE) = stats$ensembl_gene_id

pe_geneList = with(psychENCODE, 
	list(pe_SCZup = ensembl_gene_id[SCZ.t.value > 0 & SCZ.fdr < 0.05],
			pe_SCZdown = ensembl_gene_id[SCZ.t.value < 0 & SCZ.fdr < 0.05],
			pe_ASDup = ensembl_gene_id[ASD.t.value > 0 & ASD.fdr < 0.05],
			pe_ASDdown = ensembl_gene_id[ASD.t.value < 0 & ASD.fdr < 0.05],
 			pe_BDup = ensembl_gene_id[BD.t.value > 0 & BD.fdr < 0.05],
			pe_BDdown = ensembl_gene_id[BD.t.value < 0 & BD.fdr < 0.05]))

#################
## brainseq  ####
#################

## DLPFC RiboZero
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase2_qsv/rdas/dxStats_dlpfc_filtered_qSVA_noHGoldQSV_matchDLPFC.rda")

bs2_geneList = with(outGene, 
	list(bs2_SCZup = ensemblID[logFC > 0 & adj.P.Val < 0.05],
		bs2_SCZdown = ensemblID[logFC < 0 & adj.P.Val < 0.05]))


## DLPFC polyA
load("/dcl01/ajaffe/data/lab/qsva_brain/brainseq_phase1_qsv/rdas/dxStats_dlpfc_filtered_qSVA_BSP1_DLPFC.rda",verbose=TRUE)

bs1_geneList = with(outGene, 
	list(bs1_SCZup = ensemblID[logFC > 0 & adj.P.Val < 0.1],
		bs1_SCZdown = ensemblID[logFC < 0 & adj.P.Val < 0.1]))

##############################
### Sestan DS Neuron 2017? ###

###############
### combine ###
###############

## gene list ## 
geneList = c(asd_exome_geneList, asd_sfari_geneList,
	birnbaum_geneList,bs2_geneList,pe_geneList, harmonizome_geneList)

## filter for those present in spatial data
geneList_present = lapply(geneList, function(x) {
	x = x[!is.na(x)]
	x[x%in% rownames(t0_contrasts)]
})

## do enrichment
enrich_stat_list = eb0_list
for(i in seq(along=eb0_list)) {
	layer = t0_contrasts[,i] > 0 & fdrs0_contrasts[,i] < 0.1
	enrichList = mclapply(geneList_present, function(g) {
		tt = table(Set = factor(names(layer) %in% g, c(FALSE,TRUE)),
			Layer=factor(layer, c(FALSE,TRUE)))
		fisher.test(tt, alternative="greater")
	}, mc.cores=8)
	o = data.frame(OR = sapply(enrichList, "[[", "estimate"),
		Pval = sapply(enrichList, "[[", "p.value"))
	rownames(o) = gsub(".odds ratio", "", rownames(o))
	enrich_stat_list[[i]] = o
}
enrichTab = do.call("cbind", enrich_stat_list)

#  name
enrichTab$Group = rep(c("Satterstrom", "SFARI", "Birnbaum",
	"BrainSeq2", "psychENCODE",  "Harmonizome"),
	times = c(length(asd_exome_geneList), length(asd_sfari_geneList), 
			length(birnbaum_geneList), length(bs2_geneList), 
			length(pe_geneList), length(harmonizome_geneList)))
enrichTab$SetID = rownames(enrichTab) 			
enrichTab$ID = ss(rownames(enrichTab) , "_", 2)
			
enrichTab_noHarm = enrichTab[enrichTab$Group != "Harmonizome",]
pMat = enrichTab_noHarm[,grep("Pval", colnames(enrichTab_noHarm))]
colnames(pMat) = ss(colnames(pMat), "\\.")
pMat < 0.05/nrow(pMat)
pMat < 0.001

options(width=100)
signif(-log10(pMat)[c(7:9, 4:6),],3)

################
## make plots ##
################

theSeq = seq(0,12,by=0.1)					
my.col <- colorRampPalette(c("white","red"))(length(theSeq))

## overall
negPmat = -log10(pMat)
negPmat[negPmat > 12] = 12

## ASD
negPmat_ASD = negPmat[c(7:9, 4:6, 24:25),]
rownames(negPmat_ASD)[7:8] = c("DE_ASD_Up","DE_ASD_Down")

negPmat_ASD = as.matrix(negPmat_ASD)
negPmat_ASD = negPmat_ASD[,c(2:7,1)]
negPmat_ASD = negPmat_ASD[nrow(negPmat_ASD):1,]

pdf("pdf/ASD_risk_gene_heatmap.pdf",width=8)
print(levelplot(t(negPmat_ASD), aspect = "fill", at = theSeq,
	col.regions = my.col, ylab = "", xlab = "",
	scales=list(x=list(rot=90,cex=1.5), y=list(cex=1.5))))
dev.off()

enrichTab_harm = enrichTab[enrichTab$Group == "Harmonizome",]
fdrTab_harm = apply(enrichTab_harm[,grep("Pval", colnames(enrichTab_harm))], 2, p.adjust, "fdr")
colSums(fdrTab_harm < 0.05)
colSums(enrichTab_harm[,grep("Pval", colnames(enrichTab_harm))] < 1e-10)
colSums(enrichTab_harm[,grep("Pval", colnames(enrichTab_harm))] < 1e-20)
harm_top10 = apply(enrichTab_harm[,grep("Pval", colnames(enrichTab_harm))], 2, 
	function(x) gsub("Harmonizome_", "", rownames(enrichTab_harm)[order(x)[1:10]]))

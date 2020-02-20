############
library(SummarizedExperiment)
library(jaffelab)
library(readxl)
library(VariantAnnotation)
library(rtracklayer)

## load rses
load("rse_gene_hypoxia_brady_n8.Rdata")
load("rse_exon_hypoxia_brady_n8.Rdata")
load("rse_jx_hypoxia_brady_n8.Rdata")
load("rse_tx_hypoxia_brady_n8.Rdata")

# read in pheno
pd = read.csv("SraRunTable.txt",as.is=TRUE,row.names=1)
pd = pd[colnames(rse_gene),]

rse_gene$Hypoxia = ifelse(pd$source_name == "normoxia", 0 ,1)
rse_exon$Hypoxia = ifelse(pd$source_name == "normoxia", 0 ,1)
rse_jx$Hypoxia = ifelse(pd$source_name == "normoxia", 0 ,1)
rse_tx$Hypoxia = ifelse(pd$source_name == "normoxia", 0 ,1)

boxplot(mitoRate ~ Hypoxia, data = colData(rse_gene))
boxplot(concordMapRate ~ Hypoxia, data = colData(rse_gene))
boxplot(totalAssignedGene ~ Hypoxia, data = colData(rse_gene))

## drop junctions to get under file limit
rse_jx = rse_jx[rowSums(assays(rse_jx)$counts) > 1,]

## save
save(rse_gene, file="../count_data/rse_gene_hypoxia_brady_n8_annotated.Rdata")
save(rse_exon, file="../count_data/rse_exon_hypoxia_brady_n8_annotated.Rdata")
save(rse_jx, file="../count_data/rse_jx_hypoxia_brady_n8_annotated.Rdata")
save(rse_tx, file="../count_data/rse_tx_hypoxia_brady_n8_annotated.Rdata")

###########################
#### check genotypes ####
## read in VCF
vcf = readVcf("Genotypes/mergedVariants.vcf.gz")
vcf = vcf[,rse_gene$bamFile] # match up to correct BAM file
colnames(vcf) = rse_gene$SAMPLE_ID

# add rs num
snpMap = import("/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Genotyping/common_missense_SNVs_hg38.bed")
oo = findOverlaps(vcf, snpMap, type="equal")
rowData(vcf)$snpRsNum = NA
rowData(vcf)$snpRsNum[queryHits(oo)] = snpMap$name[subjectHits(oo)]
rowData(vcf)$snpRsNum[is.na(rowData(vcf)$snpRsNum)] = rownames(vcf)[is.na(rowData(vcf)$snpRsNum)]
rownames(vcf) = rowData(vcf)$snpRsNum

#########################
	######################
# subset to high-depth
vcf = vcf[info(vcf)$DP > 3*ncol(vcf) & 
          nchar(ref(vcf)) == 1 & elementNROWS(alt(vcf)) == 1,]
		 
########################################
# check snp correlation of all samples
snps = geno(vcf)$GT
snps[snps == "."] = 0
snps[snps == "0/1"] = 1
snps[snps == "1/1"] = 2
class(snps) = "numeric"
snpCor = cor(snps, use="pairwise.complete.obs")

### PLOT ####
library(pheatmap)
library(RColorBrewer)
col.pal = brewer.pal(9,"Blues")

pdf("pdfs/pheatmap_genotypes.pdf",h=10,w=10)
pheatmap(snpCor, 
		cluster_rows=TRUE, 
		cluster_cols=TRUE,
		color=col.pal)
dev.off()

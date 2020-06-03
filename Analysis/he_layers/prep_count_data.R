############
library(SummarizedExperiment)
library(jaffelab)
library(readxl)
library(VariantAnnotation)
library(rtracklayer)
library(janitor)

# read in sra and supp table info
pd = read.csv("he_SraRunTable.txt",as.is=TRUE,row.names=1)
pheno = read_excel("41593_2017_BFnn4548_MOESM254_ESM.xlsx", sheet=1)
pheno = as.data.frame( clean_names(pheno))
pheno = pheno[pheno$species == "Human",]

mm = match(pd$Sample.Name, pheno$sample_id)
pheno = pheno[mm,]
pheno$Run = pd$Run
rownames(pheno) = pheno$Run


## load rses
load("rse_gene_hypoxia_brady_n8.Rdata")
pheno = pheno[colnames(rse_gene),]
colData(rse_gene) = cbind(DataFrame(pheno), colData(rse_gene))

## save
save(rse_gene, file="rse_gene_hypoxia_brady_n8_annotated.Rdata")

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

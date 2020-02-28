#!/bin/bash
#$ -cwd
#$ -N magma_gsa
#$ -o ./logs/magma_gsa.o
#$ -e ./logs/magma_gsa.e
#$ -l mem_free=8G,h_vmem=8G

echo "**** Job starts ****"
date

model="snp-wise"
ANNO=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/MAGMA/GRCh37_ensembl_GENES_SpatialExprs.gene.loc
MAGMA=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/MAGMA/magma
BFILE=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/MAGMA/g1000_eur
PREVPATH=/dcl02/lieber/ajaffe/Nick_Clifton/magma
setcol=1
genecol=2
gs=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/MAGMA/laminar_sets.txt

mkdir -p SNP_Data
mkdir -p Results

# pgc clozuk2 schizophrenia
$MAGMA --annotate window=35,10 --snp-loc $PREVPATH/SNP_Data/clozuk_pgc2_pardinas2018.snploc --gene-loc $ANNO --out SNP_Data/clozuk_pgc2_pardinas2018_ensembl_SpatialGenes
$MAGMA --bfile $BFILE --gene-annot SNP_Data/clozuk_pgc2_pardinas2018_ensembl_SpatialGenes.genes.annot --pval $PREVPATH/SNP_Data/clozuk_pgc2_pardinas2018.meta.sumstats.txt use=SNP,P N=35802 --gene-model ${model} --out SNP_Data/clozuk_pgc2_ensembl_SpatialGenes_${model}
$MAGMA --gene-results SNP_Data/clozuk_pgc2_ensembl_SpatialGenes_snp-wise.genes.raw --set-annot $gs gene-col=${genecol} set-col=${setcol} --out Results/clozuk_pgc2_laminar

# PGC3 schizophrenia
$MAGMA --annotate window=35,10 --snp-loc $PREVPATH/SNP_Data/daner_PGC_SCZ_w3_90_0418b_INF06.snploc --gene-loc $ANNO --out SNP_Data/pgc3_scz_ensembl_SpatialGenes
$MAGMA --bfile $BFILE --gene-annot SNP_Data/pgc3_scz_ensembl_SpatialGenes.genes.annot --pval $PREVPATH/SNP_Data/daner_PGC_SCZ_w3_90_0418bN use=SNP,P ncol=Nsum --gene-model ${model} --out SNP_Data/pgc3_scz_ensembl_SpatialGenes_${model}
$MAGMA --gene-results SNP_Data/pgc3_scz_ensembl_SpatialGenes_snp-wise.genes.raw --set-annot $gs gene-col=${genecol} set-col=${setcol} --out Results/pgc3_scz_laminar

# bipolar disorder
$MAGMA --annotate window=35,10 --snp-loc $PREVPATH/SNP_Data/daner_PGC_BIP32b_mds7a_0416a_INF06.snploc --gene-loc $ANNO --out SNP_Data/daner_PGC_BIP32b_mds7a_0416a_ensembl_SpatialGenes
$MAGMA --bfile $BFILE --gene-annot SNP_Data/daner_PGC_BIP32b_mds7a_0416a_ensembl_SpatialGenes.genes.annot --pval $PREVPATH/SNP_Data/daner_PGC_BIP32b_mds7a_0416aN use=SNP,P ncol=Nsum --gene-model ${model} --out SNP_Data/daner_PGC_BIP32b_mds7a_0416a_ensembl_SpatialGenes_${model}
$MAGMA --gene-results SNP_Data/daner_PGC_BIP32b_mds7a_0416a_ensembl_SpatialGenes_snp-wise.genes.raw --set-annot $gs gene-col=${genecol} set-col=${setcol} --out Results/PGC_BIP_laminar

# depression
$MAGMA --annotate window=35,10 --snp-loc $PREVPATH/SNP_Data/MDD29_23andMe_Meta_Analysed1_CHR_POS.snploc --gene-loc $ANNO --out SNP_Data/MDD29_23andMe_Meta_Analysed1_ensembl_SpatialGenes
$MAGMA --bfile $BFILE --gene-annot SNP_Data/MDD29_23andMe_Meta_Analysed1_ensembl_SpatialGenes.genes.annot --pval $PREVPATH/SNP_Data/MDD29_23andMe_Meta_Analysed1_CHR_POS.metal use=MarkerName,P.value N=480359 --gene-model ${model} --out SNP_Data/MDD29_23andMe_ensembl_SpatialGenes_${model}
$MAGMA --gene-results SNP_Data/MDD29_23andMe_ensembl_SpatialGenes_snp-wise.genes.raw --set-annot $gs gene-col=${genecol} set-col=${setcol} --out Results/MDD29_23andMe_laminar

# autism
$MAGMA --annotate window=35,10 --snp-loc $PREVPATH/SNP_Data/iPSYCH-PGC_ASD_Nov2017.snploc --gene-loc $ANNO --out SNP_Data/iPSYCH-PGC_ASD_Nov2017_ensembl_SpatialGenes
$MAGMA --bfile $BFILE --gene-annot SNP_Data/iPSYCH-PGC_ASD_Nov2017_ensembl_SpatialGenes.genes.annot --pval $PREVPATH/SNP_Data/iPSYCH-PGC_ASD_Nov2017 use=SNP,P N=46350 --gene-model ${model} --out SNP_Data/PGC_ASD_ensembl_SpatialGenes_${model}
$MAGMA --gene-results SNP_Data/PGC_ASD_ensembl_SpatialGenes_snp-wise.genes.raw --set-annot $gs gene-col=${genecol} set-col=${setcol} --out Results/PGC_ASD_laminar

echo "**** Job ends ****"
date

### =========================================================================
### SGE variables
### -------------------------------------------------------------------------
###
#$ -l mem_free=8G,h_vmem=10G
#$ -l bluejay
#$ -m n
#$ -l h_fsize=500G
#$ -o ./logs/
#$ -e ./logs/
#$ -pe local 2
#$ -cwd
#$ -t 1-68



### =========================================================================
### Load modules
### -------------------------------------------------------------------------
###

export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2

module load python/3.6.9

conda activate ldsc3

FILELIST=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/LDSC/categories.txt
cn=$(awk "NR==$SGE_TASK_ID" $FILELIST)

### =========================================================================
### Begin code
### -------------------------------------------------------------------------
###
LDSC=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/LDSC/ldsc/ldsc.py

for sl in 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22
do


python $LDSC \
	--l2 \
	--bfile /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1/1000G_plinkfiles/1000G.mac5eur.${sl} \
	--ld-wind-cm 1 \
	--annot /dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/LDSC/annotation/${cn}.Phase1.${sl}.annot.gz \
	--out /dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/LDSC/LDScore/${cn}.Phase1.${sl} \
	--print-snps /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1/hapmap3_snps/hm.${sl}.snp \


done


# Run LDSC


### =========================================================================
### SGE variables
### -------------------------------------------------------------------------
###

#$ -l mem_free=6G
#$ -l h_vmem=8G
#$ -l h_fsize=500G
#$ -l bluejay
#$ -pe local 2
#$ -o ./logs/
#$ -e ./logs/
#$ -cwd
#$ -t 1-43


### =========================================================================
### Load modules
### -------------------------------------------------------------------------
###

module load python/3.6.9

export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2

conda activate ldsc3

FILELIST=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/LDSC/categories.txt
cn=$(awk "NR==$SGE_TASK_ID" $FILELIST)


### =========================================================================
### Adjusting for baseline
### -------------------------------------------------------------------------
###

LDSC=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/LDSC/ldsc/ldsc.py

for gwas in Autism_spectrum_disorder_latest Bipolar_disorder_latest Major_depressive_disorder_latest ADHD Agreeableness Alzheimers_disease Anorexia_nervosa Anxiety_disorder Autism_spectrum_disorder Bipolar_disorder BMI Cardioembolic_stroke Childhood_cognitive_performance Cigarettes_per_day College_attainment Conscientiousness Coronary_artery_disease Crohns_disease Depressive_symptoms Epilepsy Ever_smoked Extraversion Focal_epilepsy Generalized_epilepsy Height Intracarebral_hemorrhage IQ Ischemic_stroke Large-vessel_disease Major_depressive_disorder Neuroticism Openness PTSD Schizophrenia Small-vessel_disease Subjective_well-being Years_of_education
do


python $LDSC \
  --h2 /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/munge_sumstats/${gwas}.sumstats.gz \
  --w-ld-chr /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1/weights_hm3_no_hla/weights. \
  --ref-ld-chr /dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/LDSC/LDScore/${cn}.Phase1.,/dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1/baseline/baseline. \
  --overlap-annot \
  --frqfile-chr /dcl01/lieber/ajaffe/lab/brain-epigenomics/rdas/ldsc/Phase1/1000G_frq/1000G.mac5eur. \
  --out /dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/Layer_Guesses/LDSC/output/${cn}.${gwas}.Phase1 \
  --print-coefficients




done


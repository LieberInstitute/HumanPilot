###
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=30G,h_vmem=32G,h_fsize=100G
#$ -pe local 4
#$ -N bam_split
#$ -o logs_split_layer/bam_split.$TASK_ID.txt
#$ -e logs_split_layer/bam_split.$TASK_ID.txt
#$ -m e
#$ -t 41,46
#$ -tc 5

module load samtools
module load python/3.6.9

SUB=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X/subset-bam-1.0-x86_64-linux/subset-bam

BAMFILE=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X/sample_level_layer_map.tsv
BARFILE=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X/barcode_level_layer_map.tsv

## get bam and sample info
SAMPLE=$(awk 'BEGIN {FS="\t"} {print $1}' $BAMFILE | awk "NR==${SGE_TASK_ID}")
LAYER=$(awk 'BEGIN {FS="\t"} {print $2}' $BAMFILE | awk "NR==${SGE_TASK_ID}")
BAM=$(awk 'BEGIN {FS="\t"} {print $3}' $BAMFILE | awk "NR==${SGE_TASK_ID}")

mkdir -p /dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X/$SAMPLE/Layers/


## get out barcodes
BC_FILE=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X/$SAMPLE/Layers/${SAMPLE}_${LAYER}_barcodes.txt
awk -v a=$SAMPLE -v b=$LAYER '$2==a && $3==b {print $1}' $BARFILE > $BC_FILE


# ## subset
NEWBAM=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X/$SAMPLE/Layers/${SAMPLE}_${LAYER}.bam

$SUB --bam $BAM --cell-barcodes $BC_FILE --cores 6 --out-bam $NEWBAM

# ## index
samtools index $NEWBAM 

# ## dedup
NEWBAM_DEDUP=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X/$SAMPLE/Layers/${SAMPLE}_${LAYER}_dedup.bam

umi_tools dedup --umi-tag=UB --cell-tag=CB \
	--extract-umi-method=tag --stdin=$NEWBAM --stdout=$NEWBAM_DEDUP
	
samtools index $NEWBAM_DEDUP

## feature counts- genes
module unload conda_R

GTF=/dcl01/ajaffe/data/lab/singleCell/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf
OUTGENE=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X/$SAMPLE/Layers/${SAMPLE}_${LAYER}.genes.counts
OUTEXON=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X/$SAMPLE/Layers/${SAMPLE}_${LAYER}.exons.counts

/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/subread-1.5.0-p3-source/bin/featureCounts \
		-a $GTF -o $OUTGENE $NEWBAM_DEDUP
# exons	
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/subread-1.5.0-p3-source/bin/featureCounts  \
	-O -f -a $GTF -o $OUTEXON $NEWBAM_DEDUP


# junctions	
OUTJXN=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X/$SAMPLE/Layers/${SAMPLE}_${LAYER}.junctions.bed
OUTCOUNT=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X/$SAMPLE/Layers/${SAMPLE}_${LAYER}.junctions.count

module load python/2.7.9
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/regtools/build/regtools junctions extract -i 9 -o ${OUTJXN} ${NEWBAM_DEDUP}
/dcl01/lieber/ajaffe/Emily/RNAseq-pipeline/Software/bed_to_juncs_withCount < ${OUTJXN} > ${OUTCOUNT}

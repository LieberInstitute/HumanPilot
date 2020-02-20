#!/bin/bash
#$ -cwd
#$ -l bluejay,mf=10G,h_vmem=10G,h_fsize=100G,h_stack=256M
#$ -N download_reads
#$ -o ./logs/download.reads.$TASK_ID.txt
#$ -e ./logs/download.reads.$TASK_ID.txt
#$ -t 2-13

# https://www.ncbi.nlm.nih.gov/sra?term=SRP199498
FILELIST=/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/hafner_vglut/hafner_SraRunTable.txt
ID=$(awk "NR==$SGE_TASK_ID" $FILELIST | cut -f10)

/users/ajaffe/software/sratoolkit.2.8.1-3-centos_linux64/bin/fastq-dump \
	-O /dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/Analysis/hafner_vglut/FASTQ/ \
	--gzip --split-files $ID

#!/bin/bash

## Usage:
# sh bamtofastq.sh

## Create the logs directory
mkdir -p logs_bamtofastq

for sample in 151507 151508 151509 151510 151669 151670 151671 151672 151673 151674 151675 151676; do

    ## Internal script name
    SHORT="bamtofastq_${sample}"

    # Construct shell file
    echo "Creating script bamtofastq_${sample}"
    cat > .${SHORT}.sh <<EOF
#!/bin/bash
#$ -cwd
#$ -l bluejay,mem_free=20G,h_vmem=20G,h_fsize=100G
#$ -pe local 4
#$ -N ${SHORT}
#$ -o logs_bamtofastq/${SHORT}.txt
#$ -e logs_bamtofastq/${SHORT}.txt
#$ -m e

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: \${USER}"
echo "Job id: \${JOB_ID}"
echo "Job name: \${JOB_NAME}"
echo "Hostname: \${HOSTNAME}"
echo "Task id: \${SGE_TASK_ID}"

## List current modules for reproducibility
module list

## Edit with your job command
/dcl02/lieber/ajaffe/SpatialTranscriptomics/bamtofastq-1.2.0 --nthreads=4 /dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X/${sample}/${sample}_mRNA.bam /dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X/${sample}/fastq

echo "**** Job ends ****"
date

## This script was made using sgejobs version 0.99.1
## available from http://research.libd.org/sgejobs/


EOF

    call="qsub .${SHORT}.sh"
    echo $call
    $call
done

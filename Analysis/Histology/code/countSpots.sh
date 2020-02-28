#!/bin/bash
#$ -l mem_free=50G
#$ -l h_vmem=50G
#$ -l h_rt=24:00:00
#$ -cwd
#$ -j y
#$ -R y
#$ -t 1-12
matlab -nodisplay -nodesktop -r "addpath(genpath('/users/jcatalli/code_pipeline')); tic; try countSpots('/dcl01/lieber/ajaffe/Maddy/RNAscope/Histology/10Ximages', '/dcl02/lieber/ajaffe/SpatialTranscriptomics/HumanPilot/10X', $SGE_TASK_ID); catch e, disp(e.message); end; toc; addpath /users/jcatalli; leave"
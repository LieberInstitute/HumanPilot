#$ -l mem_free=50G,h_vmem=50G
#$ -cwd
#$ -m e
#$ -M shicks19@jhu.edu

module load conda_R/devel
python3 sce_spatialDE.py

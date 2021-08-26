#!/bin/bash
#PBS -l nodes=1:ppn=6:rhel7
#PBS -l walltime=2:00:00
#PBS -l pmem=8gb
#PBS -A mah546_c_g_bc_default #or open
#PBS -j oe
 
echo "Job started on `hostname` at `date`"

# Change working directory 
cd $PBS_O_WORKDIR

# Activate conda env
conda activate adni_project

# Set up directory
mkdir -p results
mkdir results/plots
# Copy data sets
mkdir -p data
for file in `cat ADNI_data_files.txt`
do
cp "$file" "data/"
done # Comment this for loop if manually copying the files

# Run pipeline
cd code/
echo '-----Running QC-----'
python 01_QC.py 
echo '-----Running WGCNA-----'
Rscript WGCNA.R
echo '-----Running Analysis-----'
python 02_PreliminaryAnalysis.py

echo "Job Ended at `date`"

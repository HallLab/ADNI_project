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
mkdir -p results/plots
# Copy data sets
for file in `cat ADNI_data_files.txt`
do
cp "$file" "data/"
done # Comment this for loop if manually copying the files

# Run pipeline
cd code/
echo '-----Running QC-----'
python 01_QC.py 
echo '-----Running WGCNA-----'
Rscript 02_WGCNA.R
echo '-----Running Analysis-----'
python 03_Analysis.py
echo '-----Generating Figures-----'
python 04_GenerateFigures.py
# Make sure to have imagemagick installed to merge previous plots
cd ../results/plots
montage -geometry 100% -tile 2x1 -density 150 heatmap_* Figure2.pdf
montage -geometry 100% -tile 1x2 -density 150 wgcna_power_* FigureS1.pdf

echo "Job Ended at `date`"

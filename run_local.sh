#!/bin/sh

#### CREATE DIRECTORIES
mkdir results/
mkdir results/plots

#### RUN PIPELINE
conda activate adni_project

# Run pipeline
cd code/
echo '-----Running QC-----'
python 01_QC.py 
echo '-----Running WGCNA-----'
Rscript WGCNA.R
echo '-----Running Analysis-----'
python 02_PreliminaryAnalysis.py


#!/bin/bash

#### CREATE DIRECTORIES
mkdir -p results/
mkdir -p results/plots

#### RUN PIPELINE
cd code/
echo '-----Running QC-----'
python 01_QC.py 
echo '-----Running WGCNA-----'
Rscript WGCNA.R
echo '-----Running Analysis-----'
python 02_PreliminaryAnalysis.py

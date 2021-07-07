#!/bin/sh

#### CREATE DIRECTORIES
mkdir results/

#### RUN PIPELINE
conda activate adni_project

cd code/
python 01_QC.py 
python 02_PreliminaryAnalysis.py

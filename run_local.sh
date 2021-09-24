#!/bin/bash

#### CREATE DIRECTORIES
mkdir -p results/
mkdir -p results/plots

#### RUN PIPELINE
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
montage -geometry 100% -tile 2x1 -density 150 heatmap_* Figure2.tiff
montage -geometry 100% -tile 2x1 -density 150 heatmap_* Figure2.pdf
montage -geometry 100% -tile 1x2 -density 150 wgcna_power_* FigureS1.tiff
montage -geometry 100% -tile 1x2 -density 150 wgcna_power_* FigureS1.pdf

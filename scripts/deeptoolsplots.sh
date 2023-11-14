#!/bin/bash



computeMatrix reference-point \ # choose the mode
       --referencePoint center \ # alternatives: TES, center
       -b 1500 -a 1500 \ # define the region you are interested in
       -R $1 \ #bed file
       -S $2  \ #bw file
       --skipZeros \
       -o matrix1.gz \ # to be used with plotHeatmap and plotProfile

plotHeatmap \
       -m matrix1.gz\
       -out testm1.pdf \
       --heatmapHeight 15  

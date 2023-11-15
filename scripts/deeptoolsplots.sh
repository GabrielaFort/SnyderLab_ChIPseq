#!/bin/bash



computeMatrix reference-point --referencePoint center -b 1500 -a 1500 -R $1 -S $2 --skipZeros -o matrix1.gz  # to be used with plotHeatmap and plotProfile

plotHeatmap -m matrix1.gz -out testm1.pdf --heatmapHeight 15  

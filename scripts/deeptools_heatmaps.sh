#!/bin/bash
#SBATCH -t 4:00:00 -N 1 -n 16
#SBATCH --account=snydere
#SBATCH --partition=kingspeak

# loading modules
module load bedtools
module load deeptools

# Also uses UCSC tools which are installed in the group space bin directory - must be added to path

###### Parsing command line arguments ######
arg0=$(basename "$0")

# Usage function to describe/document script
usage_info()
{
    echo "Usage: $arg0 [{-B|--bed1} condition 1 bed] [{-b|--bed2} condition 2 bed] [{-g|--genome} genome] [{-Ca|--coverage1a} condition 1 rep1 bw] [{-Cb|--coverage1b} condition 1 rep2 bw] [{-ca|--coverage2a} condition 2 rep1 bw] [-cb|--coverage2b} condition 2 rep2 bw]"
    echo -e "\nThis bash script will take two input bed files (i.e. from two different conditions)\nand four bigwig files (i.e. from two replicates of two conditions).\nIt will return heatmaps of coverage (both individual reps and merged) over the input bed files,\nindividually and merged, and will perform k-means clustering and output bed files associated with the clusters."
}


# Error function that will be called when there is no argument entered
error()
{
    echo "$arg0: $*" >&2
    exit 1
}



# Help function/message 
help_fun()
{
    usage_info
    echo
    echo "{-B|--bed1} Bed 1                  -- Condition 1 bed file"
    echo "{-b|--bed2} Bed 2                  -- Condition 2 bed file"
    echo "{-g|--genome} genome               -- Input genome that files are aligned to (mm10, mm39, hg19, hg38)"
    echo "{-Ca|--coverage1a} bw 1 rep1       -- Condition 1 Replicate 1 bigwig file"
    echo "{-Cb|--coverage1b} bw 1 rep2       -- Condition 1 Replicate 2 bigwig file"
    echo "{-ca|--coverage2a} bw 2 rep1       -- Condition 2 Replicate 1 bigwig file"
    echo "{-cb|--coverage2b} bw 2 rep2       -- Condition 2 Replicate 2 bigwig file" 
    echo "{-h|--help}                        -- Prints this help message and exits"
    exit 0
}


# Format and set input flags
# Set input args as variables (directory + genome)
# Exit and report error if any input arguments are missing 

flags()
{
    while test $# -gt 0
    do
        case "$1" in
        (-B|--bed1)
            shift
            [ $# = 0 ] && error "No condition 1 bed file specified"
            export bed1="$1"
            shift;;
        (-b|--bed2)
            shift
            [ $# = 0 ] && error "No R2 bed file specified"
            export bed2="$1"
            shift;;
        (-g|--genome)
            shift
            [ $# = 0 ] && error "No genome specified"
            export genome="$1"
            shift;;
        (-Ca|--coverage1a)
            shift
            [ $# = 0 ] && error "No condition 1, replicate 1 bigwig file specified"
            export bw1="$1"
            shift;;
        (-Cb|--coverage1b)
            shift
            [ $# = 0 ] && error "No condition 1, replicate 2 bigwig file specified"
            export bw2="$1"
            shift;;
        (-ca|--coverage2a)
            shift
            [ $# = 0 ] && error "No condition 2, replicate 1 bigwig file specified"
            export bw3="$1"
            shift;;
        (-cb|--coverage2b)
            shift
            [ $# = 0 ] && error "No condition 2, replicate 2 bigwig file specified"
            export bw4="$1"
            shift;;
        (-h|--help)
            help_fun;;
        esac
    done
}

# Run flags function on input arguments
flags "$@"


# Make sure all input arguments have been assigned to variables
if [ ! -v bed1 ]
then
  echo '-B|--bed1 is a required argument'
  exit 1
elif [ ! -v bed2 ]
then
  echo '-b|--bed2 is a required argument'
  exit 1
elif [ ! -v genome ]
then
  echo '-g|--genome is a required argument'
  exit 1
elif [ ! -v bw1 ]
then
  echo '-Ca|--coverage1a is a required argument'
  exit 1
elif [ ! -v bw2 ]
then
  echo '-Cb|--coverage1b is a required argument'
elif [ ! -v bw3 ]
then
  echo '-ca|--coverage2a is a required argument'
elif [ ! -v bw4 ]
then
  echo '-cb|--coverage2b is a required argument'
fi


# Make sure user inputs correct genome option
if [[ ! $genome =~ ^(mm10|mm39|hg19|hg38)$ ]] #use regular expressions to find either pattern
then 
  echo "Genome must be mm10, mm39, hg19, or hg38."
  exit 1
fi


# Assign chrom sizes files to variables based on input genome
if [ $genome == hg19 ]
then
  input_genome=$'/uufs/chpc.utah.edu/common/home/snydere-group1/bin/chrom_sizes_macs/hg19.chrom.sizes'
elif [ $genome == hg38 ] 
then
  input_genome=$'/uufs/chpc.utah.edu/common/home/snydere-group1/bin/chrom_sizes_macs/hg38.chrom.sizes'
elif [ $genome == mm10 ]
then
  input_genome=$'/uufs/chpc.utah.edu/common/home/snydere-group1/bin/chrom_sizes_macs/mm10.chrom.sizes'
elif [ $genome == mm39 ]
then
  input_genome=$'/uufs/chpc.utah.edu/common/home/snydere-group1/bin/chrom_sizes_macs/mm39.chrom.sizes'
fi

# Create output directory for plots
mkdir deeptools_plots

# Create merged bw files for condition 1 and 2
bigWigMerge ${bw1} ${bw2} merged.bedGraph
sort -k1,1 -k2,2n merged.bedGraph > sorted.merged.bedGraph
bedGraphToBigWig sorted.merged.bedGraph $input_genome cond1merged.bw
rm sorted.merged.bedGraph merged.bedGraph

bigWigMerge ${bw3} ${bw4} merged.bedGraph
sort -k1,1 -k2,2n merged.bedGraph > sorted.merged.bedGraph
bedGraphToBigWig sorted.merged.bedGraph $input_genome cond2merged.bw
rm sorted.merged.bedGraph merged.bedGraph

# Create merged bed files 
cat $bed1 $bed2 > union.bed
sort -k1,1 -k2,2n union.bed > sorted.union.bed
#combine overlapping peaks using bedtools
bedtools merge -i sorted.union.bed > merged.union.bed
rm union.bed sorted.union.bed


# Use deeptools to create heatmaps!

# First on two bed files and individual replicates 
computeMatrix reference-point --referencePoint center -b 1500 -a 1500 -R $bed1 $bed2 -S $bw1 $bw2 $bw3 $bw4 --missingDataAsZero --skipZeros -p 16 -o matrix.gz 
# Plot heatmap 
plotHeatmap -m matrix.gz -out deeptools_plots/peaksets_tornadoplot_reps.pdf --colorMap RdYlBu_r --legendLocation none --heatmapHeight 15 

rm matrix.gz

# Then on two bed files and combined replicates
computeMatrix reference-point --referencePoint center -b 1500 -a 1500 -R $bed1 $bed2 -S cond1merged.bw cond2merged.bw --missingDataAsZero --skipZeros -o matrix.gz
# Plot heatmap
plotHeatmap -m matrix.gz -out deeptools_plots/peaksets_tornadoplot.pdf --colorMap RdYlBu_r --legendLocation none --heatmapHeight 15

rm matrix.gz

# Now plot on merged peaksets and perform kmeans clustering
# individual replicates
computeMatrix reference-point --referencePoint center -b 1500 -a 1500 -R merged.union.bed -S $bw1 $bw2 $bw3 $bw4 --missingDataAsZero --skipZeros -p 16 -o matrix.gz
# Plot heatmap
plotHeatmap -m matrix.gz -out deeptools_plots/peakunion_tornadoplot_reps.pdf --colorMap RdYlBu_r --legendLocation none --regionsLabel Peaks --heatmapHeight 15 --outFileSortedRegions deeptools_plots/peakunion_tornadoplot_reps.bed
plotHeatmap -m matrix.gz -out deeptools_plots/peakunion_tornadoplot_reps_k2.pdf --colorMap RdYlBu_r --legendLocation none --heatmapHeight 15 --kmeans 2 --outFileSortedRegions deeptools_plots/peakunion_tornadoplot_reps_k2.bed
plotHeatmap -m matrix.gz -out deeptools_plots/peakunion_tornadoplot_reps_k4.pdf --colorMap RdYlBu_r --legendLocation none --heatmapHeight 15 --kmeans 4 --outFileSortedRegions deeptools_plots/peakunion_tornadoplot_reps_k4.bed

rm matrix.gz

# combined replicates
computeMatrix reference-point --referencePoint center -b 1500 -a 1500 -R merged.union.bed -S cond1merged.bw cond2merged.bw --missingDataAsZero --skipZeros -o matrix.gz
# Plot heatmaps
plotHeatmap -m matrix.gz -out deeptools_plots/peakunion_tornadoplot.pdf --colorMap RdYlBu_r --legendLocation none --heatmapHeight 15 --outFileSortedRegions deeptools_plots/peakunion_tornadoplot.bed
plotHeatmap -m matrix.gz -out deeptools_plots/peakunion_tornadoplot_k2.pdf --colorMap RdYlBu_r --legendLocation none --heatmapHeight 15 --kmeans 2 --outFileSortedRegions deeptools_plots/peakunion_tornadoplot_k2.bed
plotHeatmap -m matrix.gz -out deeptools_plots/peakunion_tornadoplot_k4.pdf --colorMap RdYlBu_r --legendLocation none --heatmapHeight 15 --kmeans 4 --outFileSortedRegions deeptools_plots/peakunion_tornadoplot_k4.bed

rm matrix.gz merged.union.bed cond1merged.bw cond2merged.bw







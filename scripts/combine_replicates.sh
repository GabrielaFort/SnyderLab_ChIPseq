#!/bin/bash
#SBATCH -t 6:00:00 -N 1 -n 16
#SBATCH --account=snydere
#SBATCH --partition=kingspeak

# loading modules
module load bedtools
module use $HOME/MyModules/miniconda3
# must have homer installed and in path as well! 

###### Parsing command line arguments ######
arg0=$(basename "$0")

# Usage function to describe/document script
usage_info()
{
    echo "Usage: $arg0 [{-a|--rep1} Replicate1] [{-b|--rep2} Replicate2] [{-o|--output} name] [{-g|--genome} genome]"
    echo -e "\nThis bash script will take paths to Macs2 output directories from\ntwo replicates, a reference genome, and an output directory name.\nIt will return bed files, merged bw files, and annotations and HOMER analysis on\nonly overlapping peaks between the two replicates. It will also use deeptools to\ncreate tornado plots of the two replicates side by side."
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
    echo "{-a|--rep1} Replicate1          -- Path to peak calling output directory for Rep 1"
    echo "{-b|--rep2} Replicate2          -- Path to peak calling output directory for Rep2"
    echo "{-g|--genome} genome            -- Input genome that files are aligned to (mm10, mm39, hg19, hg38)"
    echo "{-o|--output} qvalue            -- Set name of output directory" 
    echo "{-h|--help}                     -- Prints this help message and exits"
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
        (-a|--rep1)
            shift
            [ $# = 0 ] && error "No R1 bed file specified"
            export r1_path="$1"
            shift;;
        (-b|--rep2)
            shift
            [ $# = 0 ] && error "No R2 bed file specified"
            export r2_path="$1"
            shift;;
        (-g|--genome)
            shift
            [ $# = 0 ] && error "No genome specified"
            export genome="$1"
            shift;;
        (-o|--output)
            shift
            [ $# = 0 ] && error "No output directory name specified"
            export output="$1"
            shift;;
        (-h|--help)
            help_fun;;
        esac
    done
}

# Run flags function on input arguments
flags "$@"


# Make sure user inputs correct genome option
if [[ ! $genome =~ ^(mm10|mm39|hg19|hg38)$ ]] #use regular expressions to find either pattern
then 
  echo "Genome must be mm10, mm39, hg19, or hg38."
  exit 1
fi


# Make sure input directories exist
if [ ! -d $r1_path ]
then
  echo "Cannot find R1 directory OR input is not a directory"
  exit 1
elif [ ! -d $r2_path ]
then
  echo "Cannot find R2 directory OR input is not a directory"
  exit 1
fi

# Make sure all input arguments have been assigned to variables
if [ ! -v r1_path ]
then
  echo '-a|--rep1 is a required argument'
  exit 1
elif [ ! -v r2_path ]
then
  echo '-b|--rep2 is a required argument'
  exit 1
elif [ ! -v genome ]
then
  echo '-g|--genome is a required argument'
  exit 1
elif [ ! -v output ]
then
  echo '-o|--output is a required argument'
  exit 1
fi


# Find and assign files to variables
base_r1=$(basename $r1_path)
base_r2=$(basename $r2_path)
bed_r1=$(${r1_path}/${base_r1}_peaks.narrowPeak)
bed_r2=$(${r2_path}/${base_r2}_peaks.narrowPeak)
bw_r1=$(${r1_path}/${base_r1}.bw)
bw_r2=$(${r2_path}/${base_r2}.bw)
summit_r1=$(${r1_path}/${base_r1}_summits.bed)
summit_r2=$(${r2_path}/${base_r2}_summits.bed)


# Make directory with input name and navigate into it
mkdir $output
cd $output

echo -e "--------------------------------\nStarting replicate peak intersection at `date`\n------------------------------------------" > combinereps_summary.out


echo -e "Running bedtools intersect...\n" >> combinereps_summary.out
# determine largest peak file and find overlaps with other peak file
if [ wc -l ../${bed_r1} -ge wc -l ../${bed_r2} ]
then
  bedtools intersect -a ../${bed_r1} -b ../${bed_r2} -wa | uniq > ${output}_intersect.bed
  bedtools intersect -a ../${summit_r1} -b ${output}_intersect.bed -wa > ${output}_intersect_summits.bed
else
  bedtools intersect -a ../${bed_r2} -b ../${bed_r1} -wa | uniq > ${output}_intersect.bed
  bedtools intersect -a ../${summit_r2} -b ${output}_intersect.bed -wa > ${output}_intersect_summits.bed
fi

# create merged bw files
echo -e "Creating merged bw files...\n" >> combinereps_summary.out
bigWigMerge ../${bw_r1} ../${bw_r2} merged.bedGraph
bedGraphToBigWig merged.bedGraph $genome merged.bw



echo -e "Annotating peaks...\n" >> combinereps_summary.out
# Annotate peaks
annotatePeaks.pl ${output}_intersect.bed $genome -go ./GO_intersect -annStats ${output}_intersect_annotation.log > ${output}_intersect_annotation.txt

echo -e "Running homer...\n" >> combinereps_summary.out
# Run homer 
mkdir homer

# Sort bed file
sort -k1,1 -k2,2n ${output}_intersect_summits.bed | uniq | awk '{print $1,$2-50,$3+49,$4,$5}' OFS="\t" > ./homer/${output}.sorted.bed



# Now run homer to find enriched motifs
findMotifsGenome.pl ./homer/${output}.sorted.bed $genome ./homer -size 100 -mask -preparse -p 16

echo -e "Adding peak annotations to bed file...\n" >> combinereps_summary.out
# Add annotated peaks to merged bed files with python script

source $HOME/software/pkg/miniconda3/etc/profile.d/conda.sh 
# Activate conda env
conda activate chipseq

# Launch python script with appropriate command line options (will be parsed from within the script)
annotation_cleanup.py -b ${output}_intersect.bed -a ${output}_intersect_annotation.txt -g $genome

source $HOME/software/pkg/miniconda3/etc/profile.d/conda.sh
conda deactivate

num_peaks=$(wc -l ${output}_intersect.bed) 
echo -e "---------------------------------------\nJob Finished at: `date`\nNumber of intersected peaks: $num_peaks\n---------------------------------------" >> combinereps_summary.out

# Create heatmap of both replicates side by side













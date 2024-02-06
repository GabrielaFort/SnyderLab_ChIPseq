#!/bin/bash
#SBATCH -t 4:00:00 -N 1 -n 16
#SBATCH --account=snydere
#SBATCH --partition=kingspeak


# Usage function to describe/document script
usage_info()
{
    echo "Usage: $arg0 [{-d|--diffbind} Diffbind Input File] [{-f|--FDR} Optional:FDR Cutoff] [{-o|--output} name] [{-g|--genome} genome]"
    echo -e "\nThis bash script will take an input Diffbind-formatted file, an optional FDR cutoff for differential peak calling,\na reference genome, and an output directory name. It will perform Diffbind analysis and return bed files of the output and differential peaks,\nwill run HOMER on differential peaks, and will create many plots demonstrating sample comparisons."
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
    echo "{-d|--diffbind} Diffbind_input        -- Diffbind-formatted input file"
    echo "{-g|--genome} genome                  -- Input genome that files are aligned to (mm10, mm39, hg19, hg38)"
    echo "{-o|--output} name                    -- Set name of output directory"
    echo "{-f|--FDR} FDR cutoff                 -- Optional: Set FDR cutoff for differential peaks (Default=0.05)" 
    echo "{-h|--help}                           -- Prints this help message and exits"
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
        (-d|--diffbind)
            shift
            [ $# = 0 ] && error "No diffbind file specified"
            export diff_file="$1"
            shift;;
        (f-|--FDR)
            shift
            [ $# = 0 ]
            export FDR="$1"
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


# Make sure all input arguments have been assigned to variables
if [ ! -v diff_file ]
then
  echo '-d|--diffbind is a required argument'
  exit 1
elif [ ! -v FDR ]
then
  FDR="0.05"
elif [ ! -v genome ]
then
  echo '-g|--genome is a required argument'
  exit 1
elif [ ! -v output ]
then
  echo '-o|--output is a required argument'
  exit 1
fi


# Make sure user inputs correct genome option
if [[ ! $genome =~ ^(mm10|mm39|hg19|hg38)$ ]] #use regular expressions to find either pattern
then 
  echo "Genome must be mm10, mm39, hg19, or hg38."
  exit 1
fi


# Make sure diffbind file exists
if [ ! -f $diff_file ]
then
  echo "Cannot find diffbind input file!"
  exit 1
fi


# Make sure file is a csv file
if [[ ! $diff_file == *.csv ]]
then
  echo "Input file must be a csv file!"
  exit 1
fi 


# Make directory with input name and navigate into it
mkdir $output
cd $output
mv ../${diff_file} .

echo -e "-----------------------------------\nStarting Diffbind Analysis for ${diff_file}\nStarted at: `date`\n------------------------------------" > diffbind_summary.out


# Load R (version 4.1.3)
module load R/4.1.3

echo -e "\nRunning diffbind R script...\n" >> diffbind_summary.out
# Run the R script in batch, redirecting the job output to a file
diffbind.R $diff_file $FDR >> diffbind_summary.out


# Annotate output files
echo -e "Annotating diffbind output peaks...\n" >> diffbind_summary.out

# First for DEseq2 output files with all peaks 
annotatePeaks.pl diffbind_deseq_results.bed $genome > diffbind_deseq_annotate.txt
annotatePeaks.pl diffbind_edger_results.bed $genome > diffbind_edger_annotate.txt

cond1_deseq_bed=$(echo *deseq_c1_enriched.bed)
cond1_edger_bed=$(echo *edger_c1_enriched.bed)
cond2_deseq_bed=$(echo *deseq_c2_enriched.bed)
cond2_edger_bed=$(echo *edger_c2_enriched.bed)

cond1_name=$(basename $cond1_deseq_bed _deseq_c1_enriched.bed)
cond2_name=$(basename $cond2_deseq_bed _deseq_c2_enriched.bed)

# Then run on each enriched .bed file - genome is $2 - second input arg
annotatePeaks.pl $cond1_deseq_bed $genome -go ./GO_${cond1_name}_deseq -annStats ${cond1_name}_deseq_annotate.log > ${cond1_name}_deseq_annotate.txt
annotatePeaks.pl $cond1_edger_bed $genome -go ./GO_${cond1_name}_edger -annStats ${cond1_name}_edger_annotate.log > ${cond1_name}_edger_annotate.txt
annotatePeaks.pl $cond2_deseq_bed $genome -go ./GO_${cond2_name}_deseq -annStats ${cond2_name}_deseq_annotate.log > ${cond2_name}_deseq_annotate.txt
annotatePeaks.pl $cond2_edger_bed $genome -go ./GO_${cond2_name}_edger -annStats ${cond2_name}_edger_annotate.log > ${cond2_name}_edger_annotate.txt

# Run homer
echo -e "Running homer on differential peaks...\n" >> diffbind_summary.out

sort -k1,1 -k2,2n $cond1_deseq_bed | uniq > ./${cond1_name}_deseq_c1_enriched.sorted.bed
findMotifsGenome.pl ./${cond1_name}_deseq_c1_enriched.sorted.bed $genome ./${cond1_name}_deseq_homer -size given -mask -preparse

sort -k1,1 -k2,2n $cond1_edger_bed | uniq > ./${cond1_name}_edger_c1_enriched.sorted.bed
findMotifsGenome.pl ./${cond1_name}_edger_c1_enriched.sorted.bed $genome ./${cond1_name}_edger_homer -size given -mask -preparse

sort -k1,1 -k2,2n $cond2_deseq_bed | uniq > ./${cond2_name}_deseq_c2_enriched.sorted.bed
findMotifsGenome.pl ./${cond2_name}_deseq_c2_enriched.sorted.bed $genome ./${cond2_name}_deseq_homer -size given -mask -preparse

sort -k1,1 -k2,2n $cond2_edger_bed | uniq > ./${cond2_name}_edger_c2_enriched.sorted.bed
findMotifsGenome.pl ./${cond2_name}_edger_c2_enriched.sorted.bed $genome ./${cond2_name}_edger_homer -size given -mask -preparse

echo -e "Combining annotations with bed files...\n" >> diffbind_summary.out
# Make graphs and run python script to annotate all bed files...
# I want to annotate three bed files - the one with all of the results, and those enriched in each condition
# Activate conda environment
source $HOME/software/pkg/miniconda3/etc/profile.d/conda.sh 
module use $HOME/MyModules/miniconda3
conda deactivate
conda activate chipseq

annotation_cleanup.py -d diffbind_deseq_results.bed -a diffbind_deseq_annotate.txt -g $genome
annotation_cleanup.py -d diffbind_edger_results.bed -a diffbind_edger_annotate.txt -g $genome
annotation_cleanup.py -d ${cond1_deseq_bed} -a ${cond1_name}_deseq_annotate.txt -g $genome
annotation_cleanup.py -d ${cond1_edger_bed} -a ${cond1_name}_edger_annotate.txt -g $genome
annotation_cleanup.py -d ${cond2_deseq_bed} -a ${cond2_name}_deseq_annotate.txt -g $genome
annotation_cleanup.py -d ${cond2_edger_bed} -a ${cond2_name}_edger_annotate.txt -g $genome

conda activate base

rm ${cond1_name}_deseq_c1_enriched.sorted.bed ${cond2_name}_deseq_c2_enriched.sorted.bed ${cond1_name}_edger_c1_enriched.sorted.bed ${cond2_name}_edger_c2_enriched.sorted.bed

echo -e "-----------------------------------\nFinished Diffbind Analysis for ${diff_file}\nFinished at: `date`\n------------------------------------" >> diffbind_summary.out
num_peaks_c1_deseq=$(wc -l <${cond1_deseq_bed})
num_peaks_c1_edger=$(wc -l <${cond1_edger_bed})
num_peaks_c2_deseq=$(wc -l <${cond2_deseq_bed})
num_peaks_c2_edger=$(wc -l <${cond2_edger_bed})

echo -e "\nPeaks enriched in ${cond1_name} at an FDR of ${FDR}: DESEQ:${num_peaks_c1_deseq}, EDGER:${num_peaks_c1_edger}\n" >> diffbind_summary.out
echo -e "Peaks enriched in ${cond2_name} at an FDR of ${FDR}: DESEQ:${num_peaks_c2_deseq}, EDGER:${num_peaks_c2_edger}\n" >> diffbind_summary.out

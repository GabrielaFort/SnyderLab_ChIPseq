#!/bin/bash
#SBATCH -t 6:00:00 -N 1 -n 16
#SBATCH --account=snydere
#SBATCH --partition=kingspeak

# Loading macs2
module load macs/2.2.7

# must have homer installed and in path
# this script also uses bedGraphtoBigWig

###### Parsing command line arguments ######
arg0=$(basename "$0")

# Usage function to describe/document script
usage_info()
{
    echo "Usage: $arg0 [{-e|--experimental} experimental] [{-c|--control} control] [{-o|--output} name] [{-g|--genome} genome] [{-q|--qvalue} Optional:q-value]"
    echo -e "\nThis bash script will take an input experimental (chip) bam file and a control (input)\nbam file, genome, output directory name, and an optional q-value cutoff. \nIt will return Macs2 called peaks in .bed format, bigwig files for chip and control samples,\nand will run HOMER to find enriched motifs and annotate peaks."
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
    echo "{-e|--experimental} chip           -- Experimental bam file"
    echo "{-c|--control} input control       -- Control bam file"
    echo "{-o|--output} name                 -- Set name of output directory and file headers"
    echo "{-g|--genome} genome               -- Input genome that bam files are aligned to (mm10, mm39, hg19, hg38)"
    echo "{-q|--qvalue} qvalue               -- Optional: Set q-value cutoff for Macs2. Default=0.01" 
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
        (-e|--experimental)
            shift
            [ $# = 0 ] && error "No experimental sample specified"
            export chip="$1"
            shift;;
        (-c|--control)
            shift
            [ $# = 0 ] && error "No input sample specified"
            export control="$1"
            shift;;
        (-g|--genome)
            shift
            [ $# = 0 ] && error "No genome specified"
            export genome="$1"
            shift;;
        (-o|--output)
            shift
            [ $# = 0 ] && error "No output file name specified"
            export name="$1"
            shift;;
        (-q|--qvalue)
            shift
            [ $# = 0 ]
            export qvalue="$1"
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


# Make sure input files exist
if [ ! -f $chip ]
then
  echo "Cannot find experimental bam file"
  exit 1
elif [ ! -f $control ]
then
  echo "Cannot find control bam file"
  exit 1
fi

# Make sure all input arguments have been assigned to variables
if [ ! -v chip ]
then
  echo '-e|--experimental is a required argument'
  exit 1
elif [ ! -v control ]
then
  echo '-c|--control is a required argument'
  exit 1
elif [ ! -v genome ]
then
  echo '-g|--genome is a required argument'
  exit 1
elif [ ! -v name ]
then
  echo '-o|--output is a required argument'
  exit 1
elif [ ! -v qvalue ]
then
  qvalue="0.01"
fi


# Make directory with input name and navigate into it
mkdir $name
cd $name

echo "-------------------------------------------------------" > peakcalling_summary.out
echo "Starting peak calling script for $chip (experimental) and $control (control) samples" >> peakcalling_summary.out
echo "Job started at: `date`" >> peakcalling_summary.out
echo -e "-------------------------------------------------------\n" >> peakcalling_summary.out

echo -e "*******Input Parameters*******\nControl file: $control \nExperimental file: $chip \nGenome build: $genome \n qvalue: $qvalue \nName of output dir/files: $name \n" >> peakcalling_summary.out


# Assign chrom sizes files to variables based on input genome
if [ $genome == hg19 ]
then
  input_genome=$'/uufs/chpc.utah.edu/common/home/snydere-group1/bin/chrom_sizes_macs/hg19.chrom.sizes'
  gsize='hs'
elif [ $genome == hg19 ] 
then
  input_genome=$'/uufs/chpc.utah.edu/common/home/snydere-group1/bin/chrom_sizes_macs/hg38.chrom.sizes'
  gsize='hs'
elif [ $genome == mm10 ]
then
  input_genome=$'/uufs/chpc.utah.edu/common/home/snydere-group1/bin/chrom_sizes_macs/mm10.chrom.sizes'
  gsize='mm'
elif [ $genome == mm39 ]
then
  input_genome=$'/uufs/chpc.utah.edu/common/home/snydere-group1/bin/chrom_sizes_macs/mm39.chrom.sizes'
  gsize='mm'
fi


echo -e "Starting Macs2 peak calling...\n" >> peakcalling_summary.out
# Using macs2 (CHPC-installed) to call peaks
macs2 callpeak -t ../${chip} -c ../${control} --bdg -f BAMPE -g $gsize --SPMR -q $qvalue --keep-dup all -n $name

echo -e "Converting bdg files to bw files...\n" >> peakcalling_summary.out
# Creating bw file from bdg file for control and exp samples
exp_bdg=${name}_treat_pileup.bdg
ctrl_bdg=${name}_control_lambda.bdg

bedGraphToBigWig ./$exp_bdg $input_genome ./${name}.bw
bedGraphToBigWig ./$ctrl_bdg $input_genome ./${name}_control.bw

echo -e "Annotating called peaks...\n" >> peakcalling_summary.out
# Annotating called peaks using homer suite
annotatePeaks.pl ./${name}_peaks_narrowPeak $genome -go ./GO_analysis -annStats ${name}.annotation.stats.log > ${Prefix_name}_annotation.txt

echo -e "Running HOMER motif analysis...\n" >> peakcalling_summary.out
# Now running HOMER motif analysis
mkdir homer

# Use summits bed file for HOMER analysis
# Sort bed file
sort -k1,1 -k2,2n ${name}_summits.bed | uniq | awk '{print $1,$2-50,$3+49,$4,$5}' OFS="\t" > ./homer/${name}.sorted.bed

# Now run homer to find enriched motifs
findMotifsGenome.pl ./homer/${name}.sorted.bed $genome ./homer -size 100 -mask -preparse -p 16


echo -e "Updating bed file to include peak annotations:\n" >> peakcalling_summary.out
# Call python script to update narrowpeak file with annotated gene names

num_peaks=$(wc -l ${name}_peaks_narrowPeak)
echo -e "*****Summary******" >> peakcalling_summary.out
echo -e "Number of called peaks: $num_peaks \n" >> peakcalling_summary.out
echo -e "Job finished at `date`" >> peakcalling_summary.out






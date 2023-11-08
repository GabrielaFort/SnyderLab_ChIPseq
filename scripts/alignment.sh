#!/bin/bash
#SBATCH -t 24:00:00 -N 1 -n 16
#SBATCH --account=snydere
#SBATCH --partition=kingspeak


# Loading in required CHPC-installed modules

###Using Tim's UMIscripts to add SAM tags to fastq files with UMIs
module use /uufs/chpc.utah.edu/common/home/hcibcore/Modules/modulefiles
module load umiscripts

module load cutadapt
module load bowtie2
module load samtools


# Parsing command line arguments
blnk=$(echo "$arg0" | sed 's/./ /g')
arg0=$(basename "$0")

# Writing usage function 
usage_info()
{
    echo "Usage: $arg0 [{-d|--directory} path to input directory] [{-g|--genome} genome]"
    echo
    echo "This bash script will take an input directory containing fastq.gz files - each sample should have three files - R1, R2, and UMI"
    echo "Keep the same names as they are default exported from Gnomex!"
}


# Write error function that will happen when there is no argument entered
error()
{
    echo "$arg0: $*" >&2
    exit 1
}




# Write help function/documentation 
help_fun()
{
    usage_info
    echo
    echo "\n{-d|--directory} directory      -- Set path to directory with fastq files"
    echo "{-g|--genome} genome            -- Set genome to align to (mm10, mm39, hg19, hg38)"
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
        (-d|--directory)
            shift
            [ $# = 0 ] && error "No directory specified"
            export directory="$1"
            shift;;
        (-g|--genome)
            shift
            [ $# = 0 ] && error "No genome specified"
            export genome="$1"
            shift;;
        (-h|--help)
            help_fun;;
        esac
    done
}

# This runs flags function on input arguments
flags "$@"

# Make sure user inputs correct genome option
if [[ ! $genome =~ ^(mm10|mm39|hg19|hg38)$ ]] #use regular expressions to find either pattern
then 
  echo "Genome must be mm10, mm39, hg19, or hg38."
  exit 1
fi












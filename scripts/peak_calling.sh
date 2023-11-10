#!/bin/bash
#SBATCH -t 6:00:00 -N 1 -n 16
#SBATCH --account=snydere
#SBATCH --partition=kingspeak

# Loading macs2
module load macs/2.2.7

# must have homer installed and in path


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

### next need to get chrom sizes files and assign them based on input genome...





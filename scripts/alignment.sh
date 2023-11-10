#!/bin/bash
#SBATCH -t 24:00:00 -N 1 -n 16
#SBATCH --account=snydere
#SBATCH --partition=kingspeak


###### Loading in required modules ######

# Using Tim Parnell's UMIscripts and some other CHPC-installed modules
module use /uufs/chpc.utah.edu/common/home/hcibcore/Modules/modulefiles
module load umiscripts
module load cutadapt


###### Parsing command line arguments ######
arg0=$(basename "$0")

# Usage function to describe/document script
usage_info()
{
    echo "Usage: $arg0 [{-d|--directory} path to input directory] [{-g|--genome} genome]"
    echo
    echo "This bash script will take an input directory containing fastq.gz files - each sample should have three files - R1, R2, and UMI"
    echo "Keep the same names as they are default exported from Gnomex!"
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
    echo "{-d|--directory} directory      -- Set path to directory with fastq files"
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

# Run flags function on input arguments
flags "$@"


# Make sure user inputs correct genome option
if [[ ! $genome =~ ^(mm10|mm39|hg19|hg38)$ ]] #use regular expressions to find either pattern
then 
  echo "Genome must be mm10, mm39, hg19, or hg38."
  exit 1
fi


# Assign bowtie2 index paths to variables depending on input genome build...
if [ $genome == hg38 ]
then 
  input_genome=$'/uufs/chpc.utah.edu/common/home/snydere-group1/bin/bowtie_hg38_index/hg38.standard.fa'
elif [ $genome == hg19 ]
then
  input_genome=$'/uufs/chpc.utah.edu/common/home/snydere-group1/bin/bowtie_hg19_index/hg19'
elif [ $genome == mm10 ]
then
  input_genome=$'/uufs/chpc.utah.edu/common/home/snydere-group1/bin/bowtie_mm10_index/mm10.standard.bowtie2'
elif [ $genome == mm39 ]
then
  input_genome=$'/uufs/chpc.utah.edu/common/home/snydere-group1/bin/bowtie_mm39_index/GRCm39'
fi


# Make sure input directory exists
if [ ! -d $directory ]
then
  echo "Cannot find input directory."
  exit 1
fi 


# Navigate into input directory
cd $directory

echo "----------------------" > summary.out
echo "Running alignment.sh" >> summary.out
echo -e "----------------------\n" >> summary.out
echo "Run started at:" >> summary.out
date >> summary.out



# Rename files in the directory from their long gnomex names to shorter names 
for x in *gz; do mv $x ${x/_*R/_R}; done
for x in *gz; do mv $x ${x/_001/}; done


# Initiate empty array
declare -a file_names=()

# Add only sample names to empty array
for file in *gz 
do
  file_name=$(echo $file | awk -F_ '{ print $1 }')
  file_names+=($file_name)
done


# Now going to get a unique array using a dictionary
declare -A uniq_names

for filename in "${file_names[@]}"
do
  uniq_names[$filename]=0 # assigning a placeholder
done


# Iterate through unique array, find file matches, assign to R1, UMI, R3 variables, and run alignment steps on each set
for filename in "${!uniq_names[@]}"
do
  for file in *gz
  do 
    new_file_name=$(echo $file | awk -F_ '{ print $1 }')
    if [ "$new_file_name" == "$filename" ]
    then
      if [[ $file == *R1* ]]
      then 
        read_1=$file
      elif [[ $file == *R2* ]]
      then 
        UMI=$file
      elif [[ $file == *R3* ]]
      then 
        read_2=$file
      fi  
    fi
  done
  
  # Make sure all files are present for each sample - if not print error message and quit
  if [ ! -v read_1 ]
  then
    echo "An R1 file is missing for at least one sample"
    exit 1
  elif [ ! -v read_2 ]
  then
    echo "An R2 file is missing for at least one sample"
    exit 1
  elif [ ! -v UMI ]
  then
    echo "A UMI file is missing for at least one sample"
    exit 1
  fi

  echo -e "\n\n----------------------------------------" >> summary.out
  echo "Starting alignment for $filename" >> summary.out
  echo -e "----------------------------------------\n\n" >> summary.out
  echo -e "--------Appending SAM tags to fastq files with UMIs for sample $filename--------\n" >> summary.out
  

  # Use Tim's UMI scripts to add SAM tags to fastq files with UMIs
  merge_umi_fastq.pl $read_1 $read_2 $UMI &>> summary.out
  UMIfastq_1=$(basename $read_1 .fastq.gz)
  UMIfastq_2=$(basename $read_2 .fastq.gz)
  base=$(basename $read_1 _R1.fastq.gz)
  rm $read_1 $read_2
    

  ###Removing adapters from reads and quality trim
  echo -e "--------Running cutadapt on $filename to remove adaptors and quality trim reads--------\n" >> summary.out
  cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
  -o ${UMIfastq_1}.umi.cut.fastq.gz -p ${UMIfastq_2}.umi.cut.fastq.gz -q 20 -m 50 -j 16 \
  ${UMIfastq_1}.umi.fastq.gz ${UMIfastq_2}.umi.fastq.gz 
 
  rm ${UMIfastq_1}.umi.fastq.gz ${UMIfastq_2}.umi.fastq.gz

  ###Align reads using Bowtie2 (use our installed version as we need at least 2.4 and chpc's version is older
  echo -e "------------------Aligning reads and running samtools for $filename--------------------------\n"
  /uufs/chpc.utah.edu/common/home/snydere-group1/bin/bowtie2-2.4.4-linux-x86_64/bowtie2 --sam-append-comment -p 16 \
  -x ${input_genome} -1 ${UMIfastq_1}.umi.cut.fastq.gz -2 ${UMIfastq_2}.umi.cut.fastq.gz 2>> summary.out \
  | samtools fixmate -m - ${base}.bam

  samtools sort ${base}.bam -@ 32 -o ${base}.sorted.bam


  ###Using Tim's UMIscripts to discard duplicates using UMIs
  echo -e "\n------------Discarding duplicated with Tim's bam_umi_dedup.pl script for $filename-------------\n" >> summary.out
  bam_umi_dedup.pl --in ${base}.sorted.bam --distance 2500 --out ${base}.sorted.dedup.bam --cpu 12 &>> summary.out

  samtools index ${base}.sorted.dedup.bam
  

  rm ${base}.bam ${base}.sorted.bam ${base}.sorted.bam.bai ${UMI} ${UMIfastq_1}.umi.cut.fastq.gz ${UMIfastq_2}.umi.cut.fastq.gz 
  
  echo -e "\nAlignment for sample $filename finished at:" >> summary.out
  date >> summary.out
  echo >> summary.out

done


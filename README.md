# SnyderLab_ChIPseq
### Gabriela Fort and Katy Gillis

## Description
A set of scripts and tools for analysis of ChIP-sequencing datasets. Includes useful scripts for alignment
of fastq files, peak calling, peak annotation, merging replicates, and identification of differential peaks 
between samples.

## Step 1: Alignment 
### ```alignment.sh```
```
Usage: alignment.sh [{-d|--directory} path to input directory] [{-g|--genome} genome]

This bash script will take an input directory containing fastq.gz files - each sample should have three files - R1, R2, and UMI
Keep the same names as they are default exported from Gnomex!

{-d|--directory} directory      -- Set path to directory with fastq files
{-g|--genome} genome            -- Set genome to align to (mm10, mm39, hg19, hg38)
{-h|--help}                     -- Prints this help message and exits
```
This script is meant to be compatible with file names exported from Gnomex. It will take a genome build and a directory containing 
fastq files downloaded from Gnomex as input. It can handle samples from multiple experiments at once that have 'R1', 'R2', and 'R3'
in their names to designate Read 1, UMI, Read 2 fastq files (default naming from Gnomex). For each sample, this script will align
reads to the designated reference genome using [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml), deduplicate reads 
based on the provided UMI fastq files using Tim Parnell's [UMIScripts](https://github.com/HuntsmanCancerInstitute/UMIScripts/tree/master)
tool, convert sam to sorted bam files using [Samtools](http://www.htslib.org/), and will export summary files from alignment and 
deduplication.

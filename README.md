# SnyderLab_ChIPseq
### Gabriela Fort and Katy Gillis

## Description
A set of scripts and tools for analysis of ChIP-sequencing datasets. Includes useful scripts for alignment
of fastq files, peak calling, motif finding, peak annotation, merging replicates, identification of differential 
peaks between samples, and intersecting ChIP-seq and RNA-seq datasets.

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

## Step 2: Peak Calling
### ```peak_calling.sh```
```
Usage: peak_calling.sh [{-e|--experimental} experimental] [{-c|--control} control] [{-o|--output} name] [{-g|--genome} genome] [{-q|--qvalue} Optional:q-value]

This bash script will take an input experimental (chip) bam file and a control (input)
bam file, genome, output directory name, and an optional q-value cutoff. 
It will return Macs2 called peaks in .bed format, bigwig files for chip and control samples,
and will run HOMER to find enriched motifs and annotate peaks.

{-e|--experimental} chip           -- Experimental bam file
{-c|--control} input control       -- Control bam file
{-o|--output} name                 -- Set name of output directory and file headers
{-g|--genome} genome               -- Input genome that bam files are aligned to (mm10, mm39, hg19, hg38)
{-q|--qvalue} qvalue               -- Optional: Set q-value cutoff for Macs2. Default=0.01
{-h|--help}                        -- Prints this help message and exits
```


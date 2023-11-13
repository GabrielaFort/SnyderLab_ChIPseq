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

## Step 2: Peak Calling/Annotation/Motif Enrichment
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
This script will run the [Macs2](https://pypi.org/project/MACS2/) peak calling tool on two input bam files - an experimental (chip) file and a control (input) file. It will also run [HOMER motif enrichment analysis](http://homer.ucsd.edu/homer/motif/) (findMotifsGenome.pl) on called peaks and will annotate peaks to the nearest TSS and run GO enrichment analysis on annotated peaks using the HOMER suite's [annotation tool](http://homer.ucsd.edu/homer/ngs/annotation.html) and will append the annotated gene names to the output bed file from Macs2.

*Note This script relies on Python and users should follow [CHPC's instructions](https://www.chpc.utah.edu/documentation/software/python-anaconda.php) for installing miniconda into their home directories. The environment used to run these scripts must then be downloaded from the included [env.yaml](https://github.com/GabrielaFort/SnyderLab_ChIPseq/tree/main/files/chipseq.yaml) file:
```
conda env create --file chipseq.yaml --name chipseq
```

## Step 3: Combining Biological Replicates
### ```combine_replicates.sh```
```
Usage: combine_replicates.sh [{-a|--abed} Bed_R1] [{-b|--bbed} Bed_R2] [{-as|--asummit} Summit_R1] [{-bs|--bsummit} Summit_R2] [{-o|--output} name] [{-g|--genome} genome]

This bash script will take bed and summit bed files from
two replicates, a reference genome, and an output directory name.
It will return bed files, annotations and HOMER analysis on
only overlapping peaks between the two replicates.

{-a|--abed} R1_bed              -- R1 bed file
{-b|--bbed} R2_bed              -- R2 bed file
{-as|--asummit} R1_summit       -- R1 summit bed file
{-bs|--bsummit} R1_summit       -- R2 summit bed file
{-g|--genome} genome            -- Input genome that files are aligned to (mm10, mm39, hg19, hg38)
{-o|--output} qvalue            -- Set name of output directory
{-h|--help}                     -- Prints this help message and exits
```
This script uses [bedtools intersect](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) to identify overlapping peaks across two biological ChIP-seq replicates. It requires the bed files and summit files that are output from Macs2 for each biological replicate. It will return a new bed file (containing annotated gene names) of only intersected peaks, and will run HOMER motif enrichment analysis on overlapping peaks and GO analysis on overlapping annotated genes.




#!/bin/bash
#SBATCH -t 4:00:00 -N 1 -n 16
#SBATCH --account=snydere
#SBATCH --partition=kingspeak


echo "Taking input from $1"

### arg 1 is diffbind csv file
### arg 2 is genome 
### arg 3 is FDR cutoff 


mypath=$(pwd)


# Load R (version 4.1.3)
module load R/4.1.3

echo "Running diffbind R script"
# Run the R script in batch, redirecting the job output to a file
Rscript /uufs/chpc.utah.edu/common/home/snydere-group1/bin/diffbind_FDR05.R $1 $mypath > $SLURM_JOBID.out


# Annotate output files
echo "Annotating peaks"

# First for DEseq2 output files with all peaks 
annotatePeaks.pl diffbind_deseq2.txt $2 > deseq_annotate.txt
annotatePeaks.pl diffbind_edger.txt $2 > edger_annotate.txt

# Then run on each enriched .bed file - genome is $2 - second input arg
annotatePeaks.pl condition1_enriched_deseq.bed $2 -go ./GO_cond1_deseq -annStats cond1_annotate_deseq.log > cond1_annotate_deseq.txt
annotatePeaks.pl condition1_enriched_edger.bed $2 -go ./GO_cond1_edger -annStats cond1_annotate_edger.log > cond1_annotate_edger.txt

annotatePeaks.pl condition2_enriched_deseq.bed $2 -go ./GO_cond2_deseq -annStats cond2_annotate_deseq.log > cond2_annotate_deseq.txt
annotatePeaks.pl condition2_enriched_edger.bed $2 -go ./GO_cond2_edger -annStats cond2_annotate_edger.log > cond2_annotate_edger.txt


# Run homer
echo "Running homer on differential peaks"

sort -k1,1 -k2,2n condition1_enriched_deseq.bed | uniq > ./condition1_enriched_deseq.sorted.bed
findMotifsGenome.pl condition1_enriched_deseq.sorted.bed $2 ./cond1_homer_deseq -size given -mask -preparse

sort -k1,1 -k2,2n condition1_enriched_edger.bed | uniq > ./condition1_enriched_edger.sorted.bed
findMotifsGenome.pl condition1_enriched_edger.sorted.bed $2 ./cond1_homer_edger -size given -mask -preparse

sort -k1,1 -k2,2n condition2_enriched_deseq.bed | uniq > ./condition2_enriched_deseq.sorted.bed
findMotifsGenome.pl condition2_enriched_deseq.sorted.bed $2 ./cond2_homer_deseq -size given -mask -preparse

sort -k1,1 -k2,2n condition2_enriched_edger.bed | uniq > ./condition2_enriched_edger.sorted.bed
findMotifsGenome.pl condition2_enriched_edger.sorted.bed $2 ./cond2_homer_edger -size given -mask -preparse

rm condition1_enriched_deseq.bed condition1_enriched_edger.bed condition2_enriched_deseq.bed condition2_enriched_edger.bed


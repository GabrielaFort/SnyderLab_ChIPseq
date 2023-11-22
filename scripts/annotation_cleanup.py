#!/usr/bin/env python

# Must have activated chipseq conda environment 

import pandas as pd
import argparse
import sys
from argparse import RawTextHelpFormatter
import os 

# Format and document input arguments 

parser = argparse.ArgumentParser(description=f"This script cleans up the default output from homer's annotate peaks script.\nIt takes an input annotated peak file and a bed file of called peaks and\nwill append annotated gene names and other information to the bed file.", formatter_class = RawTextHelpFormatter)

parser.add_argument("-b", required = False, help='Optional: Input macs2 bed file of called peaks, if not specified, must input diffbind bed file (-d)', dest = 'bed_file')

parser.add_argument("-d", required = False, help='Optional: Input diffbind bed file, if not specified, must input macs2 bed file (-b)', dest = 'diffbind_bed_file')

parser.add_argument("-a", required = True, help='Input annotated peaks file from homer', dest = 'annotate_file')

parser.add_argument("-g", required = True, help='Input genome', dest = 'genome')

args = parser.parse_args()

if not args.annotate_file or not args.genome:
  print("Must input an annotated peak file (-a), and a genome (-g)!")
  sys.exit(1)
elif not args.bed_file and not args.diffbind_bed_file:
  print("Must input either a macs2-formatted bed file or a diffbind-formatted bed file!")
  sys.exit(1)

# assign input args to variables
annotate_file = args.annotate_file
genome = args.genome

# read in input annotation file

if args.bed_file:
  bed_file = args.bed_file
  annotate_df = pd.read_csv(annotate_file, sep="\t", header=0, names=["id","chrom","starting","ending","strnd","scre","size","annotation","det_anno","DistoTSS","promoter","Entrez","Unigene","Refseq","Ensembl","gene_name","Alias","Description","Type"])
  bed_df = pd.read_csv(bed_file, sep="\t", header=0, names=["chr","start","end","id","score","strand","signal","pval","qval","peak"])
  merge=bed_df.merge(annotate_df, how='outer', on='id')
  new=merge[["chr","start","end","id","score","strand","signal","pval","qval","peak","gene_name","Ensembl","Description"]]
  # The start, end, and score columns have been converted into floats by pandas
  # Need to change back to ints
  new=new.fillna(0)
  new=new.astype({"start":'int',"end":'int',"score":'int'})

elif args.diffbind_bed_file:
  bed_file = args.diffbind_bed_file
  annotate_df = pd.read_csv(annotate_file, sep="\t", header=0, names=["id","chr","start","ending","strnd","scre","size","annotation","det_anno","DistoTSS","promoter","Entrez","Unigene","Refseq","Ensembl","gene_name","Alias","Description","Type"])
  annotate_df["start"] = annotate_df["start"] - 1

  bed_df = pd.read_csv(bed_file, sep="\t", header=0, names=["chr","start","end","width","strand","Conc","Conc_c2","Conc_c1","Fold","p_value","FDR"])
  
  merge=bed_df.merge(annotate_df, how='outer', on=["chr","start"])
  
  
  new=merge[["chr","start","end","width","strand","Conc","Conc_c2","Conc_c1","Fold","p_value","FDR","gene_name","Ensembl","Description"]]
  new=new.fillna(0)
  new=new.astype({"start":'int',"end":'int',"width":'int'})


# save as tab separated bed file
new.to_csv(f'annotated_{bed_file}', sep='\t',index=False,header=True)














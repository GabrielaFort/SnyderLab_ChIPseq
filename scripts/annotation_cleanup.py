#!/usr/bin/env python

# Must have activated chipseq conda environment 

import pandas as pd
import argparse
import sys
from argparse import RawTextHelpFormatter
import os 

# Format and document input arguments 

parser = argparse.ArgumentParser(description=f"This script cleans up the default output from homer's annotate peaks script.\nIt takes an input annotated peak file and a bed file of called peaks and\nwill append annotated gene names and other information to the bed file.", formatter_class = RawTextHelpFormatter)

parser.add_argument("-b", required = True, help='Input bed file of called peaks', dest = 'bed_file')

parser.add_argument("-a", required = True, help='Input annotated peaks file from homer', dest = 'annotate_file')

parser.add_argument("-g", required = True, help='Input genome', dest = 'genome')

args = parser.parse_args()

if not args.bed_file or not args.annotate_file or not args.genome:
  print("Must input a bed file (-b), an annotated peak file (-a), and a genome (-g)!")
  sys.exit(1)

# assign input args to variables
bed_file = args.bed_file
annotate_file = args.annotate_file
genome = args.genome


# Read in input files as pandas dataframes
annotate_df = pd.read_csv(annotate_file, sep="\t", header=0, names=["id","chrom","starting","ending","strnd","score","size","annotation","det_anno","DistoTSS","promoter","Entrez","Unigene","Refseq","Ensembl","gene_name","Alias","Description","Type"])

bed_df = pd.read_csv(bed_file, sep="\t", header=0, names=["chr","start","end","id","score","strand","signal","pval","qval","peak"])

merge=bed_df.merge(annotate_df, how='outer', on='id')

new=merge[["chr","start","end","id","score","strand","signal","pval","qval","peak","gene_name"]]

new.to_csv(bed_file, sep='\t',index=False,header=False)



# Include mouse/human homologues in bed file next to annotated gene name
#bm = Biomart()
# note the dataset and attribute names are different
#if genome == "mm10" or genome == "mm39":
#  convert_query = bm.query(dataset='mmusculus_gene_ensembl', attributes=['external_gene_name','hsapiens_homolog_associated_gene_name'])
#  gene_homolog="hsapiens_homolog_associated_gene_name"
#elif genome == "hg19" or genome == "hg38":
#  convert_query = bm.query(dataset='hsapiens_gene_ensembl',attributes=['external_gene_name','mmusculus_homolog_associated_gene_name'])
#  gene_homolog="mmusculus_homolog_associated_gene_name"

#convert_dict = {}
#for i, row in convert_query.loc[:,["external_gene_name", gene_homolog]].iterrows():
#    if row.isna().any(): continue
#    h2m_dict[row['external_gene_name']] = row["mmusculus_homolog_associated_gene_name"]










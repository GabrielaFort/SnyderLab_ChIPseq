#!/usr/bin/env python3
# Must have activated chipseq conda environment 

import pandas as pd
import argparse
import sys
from argparse import RawTextHelpFormatter
import os
import matplotlib.pyplot as plt 
import gseapy as gp
from gseapy import Biomart


genome = sys.argv[1]
bed = sys.argv[2]

new = pd.read_csv(bed_file, sep="\t", header=0, names=["chr","start","end","id","score","strand","signal","pval","qval","peak","gene_name"])

bm=Biomart()
# Include mouse/human homologues in bed file next to annotated gene name
if genome == "mm10" or genome == "mm39":
  convert_query = bm.query(dataset='mmusculus_gene_ensembl', attributes=['external_gene_name','hsapiens_homolog_associated_gene_name'])
  gene_homolog="hsapiens_homolog_associated_gene_name"
elif genome == "hg19" or genome == "hg38":
  convert_query = bm.query(dataset='hsapiens_gene_ensembl',attributes=['external_gene_name','mmusculus_homolog_associated_gene_name'])
  gene_homolog="mmusculus_homolog_associated_gene_name"

# make dictionary of gene names and associated homologs
convert_dict = {}
for i, row in convert_query.loc[:,["external_gene_name", gene_homolog]].iterrows():
  if row.isna().any(): continue
  convert_dict[row['external_gene_name']] = row[gene_homolog]

# make new list with matching homologs and append to last column of df
homologues=[]
for gene in new["gene_name"]:
  if gene in convert_dict:
    homologues.append(gene)
  else:
    homologues.append("")
new["homolog"] = homologues



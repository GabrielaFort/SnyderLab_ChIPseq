#!/usr/bin/env python3

import pandas as pd
import sys
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
import subprocess
import argparse
from argparse import RawTextHelpFormatter
import os

###### Formatting and documenting input arguments #######

parser = argparse.ArgumentParser(description=f"This is a python script for comparison of RNAseq and ChIPseq data.\n\nIt takes a bed file from a ChIP-seq experiment with called ChIP-seq peaks and an RNAseq output file containing differential expression data.\nThis file must contain three columns - gene ID, log2FC, and p value.\nIf desired, the user may designate custom pval and log2fc cutoffs for downstream analysis. Otherwise, reasonable default values will be utilized.\nThis script will annotate called peaks based on the distance to the nearest gene TSS using HOMER.\nIt will then output many files and graphs representing overlap and intersections between significant differentially expressed genes (either up or downregulated) to called ChIP-seq peaks.", formatter_class = RawTextHelpFormatter)

parser.add_argument("-r", required = True, help='Input differential gene expression file. Required column format = gene_ID  Log2FC  p-value', dest = 'RNAseq_file')

parser.add_argument("-b", required = True, help='Supply input ChIP seq bed file', dest = 'bed_file')

parser.add_argument("-g", required = True, help='Supply genome that samples are aligned to (options: mm10, mm39, hg19, hg38', dest = 'genome')

parser.add_argument("-f", required = False, help='Optional: supply log2fc cutoff for differentially expressed gene list. Default: 0.585', dest = 'log2fc')

parser.add_argument("-p", required = False, help='Optional: supply p value cutoff for differentially expressed gene list. Default: 0.05', dest = 'pval')

args = parser.parse_args()


######### Assigning input args to variables ########

peak_file = args.bed_file
DEG_file = args.RNAseq_file
genome = args.genome

if args.pval:
  pval = args.pval
else:
  pval = 0.05

if args.log2fc:
  log2FC = args.log2fc
else:
  log2FC = 0.585


# Read in input RNAseq file
if DEG_file.endswith('.csv'):
  rna_file = pd.read_csv(DEG_file, sep=',')
elif DEG_file.endswith('.tsv'):
  rna_file = pd.read_csv(DEG_file, sep='\t')
elif DEG_file.endswith('.xls') or DEG_file.endswith('.xlsx'):
  rna_file = pd.read_excel(DEG_file)
else:
  print("Input file must be either a csv file, a tsv file, or an excel file")
  sys.exit(1)

# Export file as a tsv file without headers to parse 
rna_file.to_csv(f'new_deg_file.tsv', sep='\t', index=False, header=False)


# Annotate peaks and export results
cmd = f'annotatePeaks.pl {peak_file} {genome} > annotated_bound_genes.txt'
cmd_run = subprocess.run(cmd, shell=True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)

if cmd_run.returncode != 0:
  print(f'Peak annotation step failed')
  exit(1)


###### Functions #######
def subset_degs(deg_file, log2fc, pval):
  # This function will take a tab separated file of DEGs (gene	log2fc	pval) and will return two files
  # of up or downregulated genes based on specified log2fc and pvalue cutoffs as well as a list of sets
  # where list[0] is significantly upregulated genes and list[1] is downregulated genes

  deg_df = pd.read_csv(deg_file, sep='\t', header=None, names=['gene','log2fc','pval'])
  
  upregulated = deg_df[(deg_df['log2fc'] > float(log2fc)) & (deg_df['pval'] < float(pval))]
  upregulated.to_csv('upregulated_genes.txt',sep='\t')

  downregulated = deg_df[(deg_df['log2fc'] < -(float(log2fc))) & (deg_df['pval'] < float(pval))]
  downregulated.to_csv('downregulated_genes.txt',sep='\t')  

  up_set = set(upregulated['gene'])
  down_set = set(downregulated['gene'])
  
  set_list = []
  set_list.append(up_set)
  set_list.append(down_set)
  
  return set_list


def set_intersection(set_list, anno_chip):
  # This function will take a list of sets of DEGs where list[0] is upregulated genes and list[1] is
  # downregulated genes. It will also take a file containing annotated chip seq peaks where the gene
  # names are in column[15] (default output from HOMER's peak annotation). 
  # The function will return a list of sets where list[0] is upregulated bound genes and list[1] is
  # downregulated bound genes
  chip_df = pd.read_csv(anno_chip, sep='\t', header = 0)
  peak_set = set(chip_df['Gene Name'])
  
  intersect_list = []
  
  for DE_set in set_list:
    intersect = DE_set.intersection(peak_set)
    intersect_list.append(intersect)

  return intersect_list


def venn_diagram(set1, set2, set3, name1, name2, name3):
  # This function will take three sets and three names and will create a three-way venn diagram
  # showing the intersection between all three sets with the appropriate name labels
  venn3([set1,set2,set3], (name1, name2, name3))
  plt.savefig('venn_diagram.pdf')


def bar_graph(x_list, y_list):
  # This function will take two lists, one with X values and one with Y values, and will make a bar graph
  fig,ax = plt.subplots(figsize = (8,8))
  ax.bar(x_list, y_list, color = 'maroon', edgecolor = 'black')
  ax.set_ylabel("% DEGs bound by TF")
  ax.set_title("RNA-seq Differential Genes vs ChIP-seq Peaks")
  fig.savefig('bar_graph.pdf')


# Run function to subset input DEG file
deg_set_list = subset_degs('new_deg_file.tsv', log2FC, pval)

# Intersect annotated peaks and DEGs
intersect_list = set_intersection(deg_set_list, 'annotated_bound_genes.txt')

# create output text files containing up or downregulated genes bound by TF
with open('upregulated_bound_genes.txt', 'w') as up_write:
  for bound_gene in intersect_list[0]:
    up_write.write(f'{bound_gene}\n')
  
with open('downregulated_bound_genes.txt','w') as down_write:
  for bound_gene in intersect_list[1]:
    down_write.write(f'{bound_gene}\n')


with open('intersect_summary.txt', 'w') as summary_write:
  summary_write.write(f'_____________________________________________________________________________________________________________________\n')
  summary_write.write(f'''Number of upregulated genes: {len(deg_set_list[0])}\nNumber of downregulated genes: {len(deg_set_list[1])}\n\n
Number of upregulated genes bound by TF: {len(intersect_list[0])}\nNumber of downregulated genes bound by TF: {len(intersect_list[1])}\n\n
Fraction upregulated genes bound by TF: {len(intersect_list[0])/len(deg_set_list[0])}\nFraction downregulated genes bound by TF: {len(intersect_list[1])/len(deg_set_list[1])}''')
  summary_write.write(f'\n_____________________________________________________________________________________________________________________\n')



# create a few different output graphs showing overlap between degs and TF-bound genes!

# First making a venn diagram showing intersections
upregulated_genes_set = deg_set_list[0]
downregulated_genes_set = deg_set_list[1]
chip_df = pd.read_csv('annotated_bound_genes.txt', sep='\t', header = None)
bound_genes_set = set(chip_df[15])

venn_diagram(upregulated_genes_set, downregulated_genes_set, bound_genes_set, 'Up_Genes','Down_Genes','Bound_Genes')



# Now making a bar graph showing intersections
percent_up_bound = len(intersect_list[0])/len(deg_set_list[0])
percent_down_bound = len(intersect_list[1])/len(deg_set_list[1])

y_list = [percent_up_bound,percent_down_bound]
x_list = ['Up Genes','Down Genes']

bar_graph(x_list, y_list)


os.remove('new_deg_file.tsv')
os.remove('annotated_bound_genes.txt')


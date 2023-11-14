#!/usr/bin/env python

# Must have activated chipseq conda environment 

import pandas as pd
import argparse
import sys
from argparse import RawTextHelpFormatter
import os 

# Format and document input arguments 

parser = argparse.ArgumentParser(description=f"This script plots heatmaps and signal plots for either individual samples or replicates.\nIt takes either one or two bed/bw files\nwill append annotated gene names and other information to the bed file.", formatter_class = RawTextHelpFormatter)

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


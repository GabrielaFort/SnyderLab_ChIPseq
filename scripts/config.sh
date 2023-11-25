#!/bin/bash

# Get and install miniconda3 into home directory
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash ./Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/software/pkg/miniconda3 -s

# Create a user environment module
mkdir -p $HOME/MyModules/miniconda3
cp /uufs/chpc.utah.edu/sys/installdir/python/modules/miniconda3/latest.lua $HOME/MyModules/miniconda3

# Load chipseq conda environment for chipseq analysis
module use $HOME/MyModules
source $HOME/software/pkg/miniconda3/etc/profile.d/conda.sh

conda env create --file ../files/chipseq.yaml --name chipseq

conda init

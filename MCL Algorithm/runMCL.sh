#!/bin/bash

# Input file
graph="nodes.txt"

# Activate virtual environment
conda activate

# call MCL algorithm
mcl ${graph} --abc -I 2 -o kroganMCL_2.txt &> kroganOutput_2.txt

# close virtual environment
conda deactivate

#!/bin/bash

# Input file
graph="nodes.txt"

conda activate

for i in 1.5
do
mcl ${graph} --abc -I ${i} -o kroganMCL_${i}.txt &> kroganOutput_${i}.txt
done

conda deactivate


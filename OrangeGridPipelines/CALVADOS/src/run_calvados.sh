#!/bin/bash

trial=$1
protein=$2

# Change to the correct directory
cd /home/$(whoami)/CALVADOS/src/
eval "$(/home/$(whoami)/miniconda3/bin/conda shell.bash hook)"
conda activate CALVADOS3

export RAY_DISABLE_MEMORY_MONITOR=1

echo "Trial: $trial, Protein: $protein"
echo "CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"
echo "Python: $(which python)"

# Run Python job
python single_IDR.py "$trial" "$protein"

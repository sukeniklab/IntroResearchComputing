#!/bin/bash

PROTEIN_NAME=$1
MASTER_FASTA="$HOME/batch_nardini/fastas/homo_sapiens_short.fasta"
OUTPUT_DIR="$HOME/batch_nardini/output_csvs"


mkdir -p "$OUTPUT_DIR"

# activate conda (two-step pattern needed since HTCondor jobs run non-interactively)
eval "$(/home/$(whoami)/miniconda3/bin/conda shell.bash hook)"
conda activate nardini_features  

python "$HOME/batch_nardini/calc_features_single.py" \
    "$PROTEIN_NAME" \
    "$MASTER_FASTA" \
    "$OUTPUT_DIR" 

echo "Done: $PROTEIN_NAME"

#!/bin/bash

FASTA_FILE=$1
BASENAME=$(basename "$FASTA_FILE" .fasta)
SPECIES_NAME=$(echo "$BASENAME" | sed -E 's/_chunk[0-9]+$//')
OUTPUT_DIR="$HOME/esm_embed/output_pt/${SPECIES_NAME}"

echo "Starting: $SPECIES_NAME"
echo "FASTA: $FASTA_FILE"
echo "Output: $OUTPUT_DIR"

# activate conda
eval "$(/home/$(whoami)/miniconda3/bin/conda shell.bash hook)"
conda activate esm2


mkdir -p "$OUTPUT_DIR"

python "$HOME/esm_embed/esm/scripts/extract.py" esm2_t33_650M_UR50D \
    "$FASTA_FILE" \
    "$OUTPUT_DIR" \
    --repr_layers $(seq 0 33) \
    --include mean per_tok \
    --toks_per_batch 4096 \
    --truncation_seq_length 3000

echo "Done: $SPECIES_NAME"


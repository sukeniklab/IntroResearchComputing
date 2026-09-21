#!/bin/bash

PDB_FILE=$1
SEQUENCE=$(basename "$PDB_FILE" .pdb)
WORK_DIR="/home/$(whoami)/CALVADOS/init_AF_struc/${SEQUENCE}"
OUTPUT_DIR="/home/$(whoami)/CALVADOS/src/starting"

echo "Starting: $SEQUENCE"
echo "PDB: $PDB_FILE"
echo "Working dir: $WORK_DIR"

# activate conda
eval "$(/home/$(whoami)/miniconda3/bin/conda shell.bash hook)"
conda activate modeller_env

# Each job gets its own working directory to avoid file collisions
mkdir -p "$WORK_DIR"
cp "$PDB_FILE" "$WORK_DIR/"
cd "$WORK_DIR"

python run_modeller.py "SEQUENCE"

# Rename modeller outputs from target.B9999000X.pdb -> {sequence}_trial{n}.pdb
trial=1
for f in $(ls target.B*.pdb 2>/dev/null | sort); do
    mv "$f" "${OUTPUT_DIR}/${SEQUENCE}_trial${trial}.pdb"
    echo "Wrote: ${SEQUENCE}_trial${trial}.pdb"
    ((trial++))
done
mv "${SEQUENCE}-seq.txt" "${OUTPUT_DIR}/${SEQUENCE}-seq.txt"

echo "Done: $SEQUENCE"
#!/usr/bin/env python3
# split_fastas.py
import sys, os
from pathlib import Path
from Bio import SeqIO

SRC_DIR = Path(".")
OUT_DIR = Path("to_run")
CHUNK_SIZE = 50  ## smaller = more jobs submitted at once 

for fasta_path in SRC_DIR.glob("*.fasta"):  ## for each fasta file
    species = fasta_path.stem
    records = list(SeqIO.parse(fasta_path, "fasta"))

    seen_ids = {}
    header_map_path = OUT_DIR / f"{species}_header_map.tsv"
    with open(header_map_path, "w") as map_f:
        for record in records:
            ## clean the sequence headers 
            ## trims excess from uniprot
            safe_id = record.id.translate(str.maketrans("/\\:", "___"))
            if safe_id in seen_ids:
                ## duplicate check 
                seen_ids[safe_id] += 1
                safe_id = f"{safe_id}_dup{seen_ids[safe_id]}"
            else:
                seen_ids[safe_id] = 0
            map_f.write(f"{safe_id}\t{record.description}\n")
            record.id = safe_id
            record.name = safe_id
            record.description = safe_id  

    ## write the chunk fastas
    for i in range(0, len(records), CHUNK_SIZE): 
        chunk = records[i:i+CHUNK_SIZE]
        chunk_idx = i // CHUNK_SIZE
        chunk_path = OUT_DIR/ f"{species}_chunk{chunk_idx:04d}.fasta"
        SeqIO.write(chunk, chunk_path, "fasta")


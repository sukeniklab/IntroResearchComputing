#!/usr/bin/env python3


import argparse


def extract_accession(header):
    parts = header.split('|')
    if len(parts) >= 3:
        return parts[1].strip()
    return header.strip()


def main():
    parser = argparse.ArgumentParser(
        description="Extract UniProt accessions from a FASTA into a flat, one-per-line text file."
    )
    parser.add_argument("master_fasta", help="Path to the full FASTA")
    parser.add_argument("output_txt", help="Path to write the flat list of accessions to")
    args = parser.parse_args()

    accessions = []
    with open(args.master_fasta) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                accessions.append(extract_accession(line[1:]))

    with open(args.output_txt, "w") as f:
        for acc in accessions:
            f.write(acc + "\n")

    n_unique = len(set(accessions))
    print(f"Found {len(accessions)} sequences in {args.master_fasta}")
    if n_unique != len(accessions):
        print(f"  WARNING: only {n_unique} unique accessions -- {len(accessions) - n_unique} duplicate(s) found")
    print(f"Wrote list to {args.output_txt}")


if __name__ == "__main__":
    main()

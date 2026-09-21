#!/usr/bin/env python3

import argparse
import glob
import os
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Combine per-sequence CSVs into one final CSV.")
    parser.add_argument("csv_dir", help="Directory of one-row CSVs, one per sequence")
    parser.add_argument("output_csv", help="Path to write the final combined CSV")
    args = parser.parse_args()

    csv_files = sorted(glob.glob(os.path.join(args.csv_dir, "*.csv")))
    print(f"Found {len(csv_files)} completed sequence CSVs in {args.csv_dir}")

    if not csv_files:
        raise ValueError("No CSVs found -- check that all condor jobs finished.")

    all_rows = [pd.read_csv(f) for f in csv_files]
    final_df = pd.concat(all_rows, ignore_index=True)
    final_df.to_csv(args.output_csv, index=False)
    print(f"Wrote {len(final_df)} total sequences to {args.output_csv}")


if __name__ == "__main__":
    main()

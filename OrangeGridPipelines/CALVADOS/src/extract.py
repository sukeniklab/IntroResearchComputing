#!/usr/bin/env python3
"""
For each fusion protein (triplicate), find the trajectory frame with the
smallest radius of gyration and write it out as a PDB.

Assumes pre-computed Rg_traj.npy files. Falls back to computing Rg on the
fly with MDAnalysis if the .npy is absent.

Usage:
    python extract_compact_fusion.py \
        --csv  /path/to/multidomain_metapredictv3_Dec17.csv \
        --base /home/jkniblo/Fusion/Data/fusions \
        --out  /path/to/output_dir \
        [--trials 3]
"""

import argparse
import ast
import csv
import os
import warnings

import numpy as np

try:
    import MDAnalysis as mda
    from MDAnalysis.analysis.rms import RMSD
except ImportError:
    raise ImportError("MDAnalysis is required: pip install MDAnalysis")


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def build_paths(base: str, fusion_name: str, trial: int):
    """Return (sim_dir, pdb_path, dcd_path, rg_npy_path) for one run."""
    sim_dir = os.path.join(
        base,
        f"CALVADOS3COM_2.0_MD_gpu_trial{trial}_{fusion_name}",
        fusion_name,
        "0",
    )
    pdb  = os.path.join(sim_dir, f"{fusion_name}.pdb")
    dcd  = os.path.join(sim_dir, f"{fusion_name}.dcd")
    npy  = os.path.join(sim_dir, "Rg_traj.npy")
    return sim_dir, pdb, dcd, npy


def get_rg_array(pdb: str, dcd: str, npy: str) -> np.ndarray:
    """Load pre-computed Rg or compute it frame-by-frame."""
    if os.path.exists(npy):
        return np.load(npy)

    warnings.warn(f"Rg_traj.npy not found, computing from trajectory: {dcd}")
    u = mda.Universe(pdb, dcd)
    rg_vals = []
    for _ in u.trajectory:
        rg_vals.append(u.atoms.radius_of_gyration())
    return np.array(rg_vals)


def extract_min_rg_frame(pdb: str, dcd: str, npy: str, out_pdb: str):
    """Write the most compact frame to out_pdb."""
    rg_vals = get_rg_array(pdb, dcd, npy)
    min_idx = int(np.argmin(rg_vals))
    min_rg  = rg_vals[min_idx]

    u = mda.Universe(pdb, dcd)
    u.trajectory[min_idx]
    u.atoms.write(out_pdb)

    return min_idx, min_rg


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv",     required=False, help="Input TSV with fusion_name and fd_boundaries columns")
    parser.add_argument("--protein", required=False, help="Single fusion name to run instead of --csv (e.g. RTL8C_FL)")
    parser.add_argument("--base",   required=True, help="Base directory containing CALVADOS run folders")
    parser.add_argument("--out",    required=True, help="Output directory for extracted PDBs")
    parser.add_argument("--trials", type=int, default=3, help="Number of replicates (default: 3)")
    args = parser.parse_args()

    if not args.csv and not args.protein:
        parser.error("must provide either --csv or --protein")

    os.makedirs(args.out, exist_ok=True)

    # Read fusion names from CSV/TSV, or use the single protein given directly
    if args.protein:
        fusion_names = [args.protein.strip()]
    else:
        fusion_names = []
        with open(args.csv, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                fusion_names.append(row["fusion_name"].strip())

    print(f"Found {len(fusion_names)} fusions × {args.trials} trials = "
          f"{len(fusion_names) * args.trials} total runs\n")

    missing  = []
    errors   = []
    summary  = []   # (fusion, trial, min_idx, min_rg, out_pdb)

    for fusion in fusion_names:
        for trial in range(1, args.trials + 1):
            sim_dir, pdb, dcd, npy = build_paths(args.base, fusion, trial)

            # Check required files exist
            if not os.path.exists(pdb) or not os.path.exists(dcd):
                missing.append(f"{fusion} trial{trial}  →  {sim_dir}")
                continue

            out_pdb = os.path.join(args.out, f"{fusion}_trial{trial}_COM_ini.pdb")

            try:
                min_idx, min_rg = extract_min_rg_frame(pdb, dcd, npy, out_pdb)
                print(f"[OK]  {fusion:40s}  trial{trial}  "
                      f"frame={min_idx:5d}  Rg={min_rg:.3f} Å  →  {out_pdb}")
                summary.append((fusion, trial, min_idx, round(float(min_rg), 4), out_pdb))
            except Exception as exc:
                errors.append(f"{fusion} trial{trial}: {exc}")
                print(f"[ERR] {fusion} trial{trial}: {exc}")

    # Write a summary TSV
    summary_path = os.path.join(args.out, "minRg_summary.tsv")
    with open(summary_path, "w") as fh:
        fh.write("fusion_name\ttrial\tmin_frame\tmin_Rg_A\tout_pdb\n")
        for row in summary:
            fh.write("\t".join(str(x) for x in row) + "\n")

    print(f"\nDone.  Summary written to {summary_path}")

    if missing:
        print(f"\n⚠  {len(missing)} run(s) skipped (files not found):")
        for m in missing:
            print(f"   {m}")

    if errors:
        print(f"\n✗  {len(errors)} run(s) failed:")
        for e in errors:
            print(f"   {e}")


if __name__ == "__main__":
    main()

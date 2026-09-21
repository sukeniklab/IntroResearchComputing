#!/usr/bin/env python3
"""
Unified IDP Analysis: Nardini+ features, CIDER parameters, and SPARROW parameters
Adapted from nardini_plus_final.ipynb for use as a standalone, per-chunk script.
Same calculation logic as the notebook — only the I/O layer changed (argparse
instead of hardcoded USER INPUTS, and no display() calls).
"""

import argparse
import numpy as np
import pandas as pd
import re
from localcider.sequenceParameters import SequenceParameters as CIDER_SP
from sparrow import Protein

def main():
    parser = argparse.ArgumentParser(description="Compute Nardini+/CIDER/SPARROW features for one fasta chunk")
    parser.add_argument("input_fasta", help="Path to input fasta chunk")
    parser.add_argument("output_csv", help="Path to output CSV")
    parser.add_argument("--calculate_zscores", action="store_true", default=False)
    parser.add_argument("--reference_species", default="Homo sapiens")
    args = parser.parse_args()

    input_fasta_file = args.input_fasta
    calculate_zscores = args.calculate_zscores
    reference_species = args.reference_species

    # ========================================================================
    # LOAD REFERENCE DATA FOR Z-SCORES
    # ========================================================================
    if calculate_zscores:
        print("="*70)
        print("Z-score calculation enabled.")
        print(f"Loading reference data for: {reference_species}")
        print("="*70)

        sheetID = '1yxt0R1G0gdI2bGpjYXgk_h7-1qA6EY6J'
        worksheetName = reference_species.split(" ")[1]
        currurl = f'https://docs.google.com/spreadsheets/d/{sheetID}/gviz/tq?tqx=out:csv&sheet={worksheetName}'

        try:
            speciesdf = pd.read_csv(currurl)
            meanvals_species = speciesdf['Mean'].values
            stdvals_species = speciesdf['Std'].values
            mycompfeats_all = speciesdf['Feature'].tolist()
            print(f"Successfully loaded {len(mycompfeats_all)} Nardini features.\n")
        except Exception as e:
            raise ValueError(f"Failed to load reference data: {e}")
    else:
        print("="*70)
        print("Z-score calculation disabled.")
        print("="*70 + "\n")

        aas = 'ACDEFGHIKLMNPQRSTVY'
        mycompfeats_all = ['fracA', 'fracC', 'fracD', 'fracE', 'fracF', 'fracG', 'fracH', 'fracI',
                            'fracK', 'fracL', 'fracM', 'fracN', 'fracP', 'fracQ', 'fracR', 'fracS',
                            'fracT', 'fracV', 'fracW', 'fracY', 'fracpos', 'fracneg', 'fracpol',
                            'fracali', 'fracaro', 'fracRtoK', 'fracEtoD', 'fracexp', 'fcr', 'ncpr',
                            'mhydro', 'dispro', 'isopoi', 'ppii']
        for aa in aas:
            mycompfeats_all.append(f'patch{aa}')
        mycompfeats_all.append('patchRG')

    # ========================================================================
    # LOAD FASTA FILE
    # ========================================================================
    print(f"Loading FASTA file: {input_fasta_file}")
    with open(input_fasta_file, 'r') as myfile:
        Lines = myfile.readlines()

    subseqs = []
    subnames = []
    thisseq = ''

    for line in Lines:
        cleanline = line.strip()
        if cleanline.startswith('>'):
            subnames.append(cleanline[1:])
            if thisseq != '':
                subseqs.append(thisseq.upper())
                thisseq = ''
        else:
            thisseq += cleanline

    subseqs.append(thisseq.upper())
    print(f"Loaded {len(subseqs)} sequences\n")

    # ========================================================================
    # CALCULATE NARDINI+ COMPOSITIONAL FEATURES
    # ========================================================================
    numInt = 2
    minBlockLen = 4
    aas = 'ACDEFGHIKLMNPQRSTVY'

    fracA, fracC, fracD, fracE, fracF, fracG, fracH, fracI, fracK, fracL = [], [], [], [], [], [], [], [], [], []
    fracM, fracN, fracP, fracQ, fracR, fracS, fracT, fracV, fracW, fracY = [], [], [], [], [], [], [], [], [], []
    fracpos, fracneg, fracpol, fracali, fracaro = [], [], [], [], []
    fracRtoK, fracEtoD, fracexp, fcr, ncpr = [], [], [], [], []
    mhydro, dispro, isopoi, ppii = [], [], [], []
    fracpatch = [[] for _ in range(len(aas))]
    rgpatch = []

    print("Calculating Nardini+ compositional features...")

    for currseq in subseqs:
        if len(currseq) >= 1 and not any(x in currseq for x in "XUZJBO"):
            SeqOb = CIDER_SP(currseq)
            slen = SeqOb.get_length()
            aafrac = SeqOb.get_amino_acid_fractions()

            fracexp.append(SeqOb.get_fraction_expanding())
            fcr.append(SeqOb.get_FCR())
            ncpr.append(SeqOb.get_NCPR())
            mhydro.append(SeqOb.get_mean_hydropathy())
            dispro.append(SeqOb.get_fraction_disorder_promoting())
            isopoi.append(SeqOb.get_isoelectric_point())
            ppii.append(SeqOb.get_PPII_propensity(mode='hilser'))

            fracA.append(aafrac['A']); fracC.append(aafrac['C']); fracD.append(aafrac['D'])
            fracE.append(aafrac['E']); fracF.append(aafrac['F']); fracG.append(aafrac['G'])
            fracH.append(aafrac['H']); fracI.append(aafrac['I']); fracK.append(aafrac['K'])
            fracL.append(aafrac['L']); fracM.append(aafrac['M']); fracN.append(aafrac['N'])
            fracP.append(aafrac['P']); fracQ.append(aafrac['Q']); fracR.append(aafrac['R'])
            fracS.append(aafrac['S']); fracT.append(aafrac['T']); fracV.append(aafrac['V'])
            fracW.append(aafrac['W']); fracY.append(aafrac['Y'])

            fracpos.append(aafrac['K'] + aafrac['R'])
            fracneg.append(aafrac['D'] + aafrac['E'])
            fracpol.append(aafrac['Q'] + aafrac['N'] + aafrac['S'] + aafrac['T'] + aafrac['G'] + aafrac['C'] + aafrac['H'])
            fracali.append(aafrac['A'] + aafrac['L'] + aafrac['M'] + aafrac['I'] + aafrac['V'])
            fracaro.append(aafrac['F'] + aafrac['W'] + aafrac['Y'])

            fracRtoK.append(np.log10(((slen * aafrac['R']) + 1) / ((slen * aafrac['K']) + 1)))
            fracEtoD.append(np.log10(((slen * aafrac['E']) + 1) / ((slen * aafrac['D']) + 1)))

            for idx, aa in enumerate(aas):
                justKs = '0' * len(currseq)
                pos = [i for i, ltr in enumerate(currseq) if ltr == aa]

                pos2 = pos.copy()
                for p in range(len(pos) - 1):
                    tdi = pos[p + 1] - pos[p]
                    if 1 < tdi <= numInt + 1:
                        pos2.extend(range(pos[p] + 1, pos[p + 1]))

                justKs = list(justKs)
                for p in pos2:
                    justKs[p] = '1'
                justKs = ''.join(justKs)

                the_ones = re.findall(r"1+", justKs)
                idx_ones = [[m.start(0), m.end(0)] for m in re.finditer(r"1+", justKs)]

                patchescombined = ''
                for count, o in enumerate(the_ones):
                    myrange = idx_ones[count]
                    subseq = currseq[myrange[0]:myrange[1]]
                    pos3 = [i for i, ltr in enumerate(subseq) if ltr == aa]
                    if len(pos3) >= minBlockLen:
                        patchescombined += subseq

                fracpatch[idx].append(len(patchescombined) / len(currseq))

            justKs = '0' * len(currseq)
            pos = [i for i, ltr in enumerate(currseq) if ltr in 'RG']

            pos2 = pos.copy()
            for p in range(len(pos) - 1):
                tdi = pos[p + 1] - pos[p]
                if 1 < tdi <= numInt + 1:
                    pos2.extend(range(pos[p] + 1, pos[p + 1]))

            justKs = list(justKs)
            for p in pos2:
                justKs[p] = '1'
            justKs = ''.join(justKs)

            the_ones = re.findall(r"1+", justKs)
            idx_ones = [[m.start(0), m.end(0)] for m in re.finditer(r"1+", justKs)]

            patchescombined = ''
            for count, o in enumerate(the_ones):
                myrange = idx_ones[count]
                subseq = currseq[myrange[0]:myrange[1]]
                if subseq.count('RG') >= 2:
                    patchescombined += subseq

            rgpatch.append(len(patchescombined) / len(currseq))

        else:
            fracexp.append(np.nan); fcr.append(np.nan); ncpr.append(np.nan)
            mhydro.append(np.nan); dispro.append(np.nan); isopoi.append(np.nan); ppii.append(np.nan)
            fracA.append(np.nan); fracC.append(np.nan); fracD.append(np.nan); fracE.append(np.nan)
            fracF.append(np.nan); fracG.append(np.nan); fracH.append(np.nan); fracI.append(np.nan)
            fracK.append(np.nan); fracL.append(np.nan); fracM.append(np.nan); fracN.append(np.nan)
            fracP.append(np.nan); fracQ.append(np.nan); fracR.append(np.nan); fracS.append(np.nan)
            fracT.append(np.nan); fracV.append(np.nan); fracW.append(np.nan); fracY.append(np.nan)
            fracpos.append(np.nan); fracneg.append(np.nan); fracpol.append(np.nan)
            fracali.append(np.nan); fracaro.append(np.nan)
            fracRtoK.append(np.nan); fracEtoD.append(np.nan)
            for idx in range(len(aas)):
                fracpatch[idx].append(np.nan)
            rgpatch.append(np.nan)

    compfeatvals_all = [fracA, fracC, fracD, fracE, fracF, fracG, fracH, fracI, fracK, fracL,
                        fracM, fracN, fracP, fracQ, fracR, fracS, fracT, fracV, fracW, fracY,
                        fracpos, fracneg, fracpol, fracali, fracaro, fracRtoK, fracEtoD,
                        fracexp, fcr, ncpr, mhydro, dispro, isopoi, ppii]

    for a in fracpatch:
        compfeatvals_all.append(a)
    compfeatvals_all.append(rgpatch)

    print("Nardini+ features calculated.\n")

    # ========================================================================
    # CREATE NARDINI DATAFRAME WITH RAW VALUES AND Z-SCORES
    # ========================================================================
    num_seqs = len(subseqs)
    num_feats = len(mycompfeats_all)
    subrawcomp = np.zeros((num_seqs, num_feats))

    for s in range(num_seqs):
        for f in range(num_feats):
            subrawcomp[s, f] = compfeatvals_all[f][s]

    raw_df = pd.DataFrame(data=subrawcomp, columns=mycompfeats_all)

    if calculate_zscores:
        print("Calculating z-scores for Nardini features...")
        subzveccomp = np.zeros((num_seqs, num_feats))
        min_std_threshold = 1e-10
        zero_std_features = [mycompfeats_all[f] for f in range(num_feats) if stdvals_species[f] < min_std_threshold]

        if zero_std_features:
            print(f"Warning: {len(zero_std_features)} feature(s) have zero/near-zero std dev")

        for s in range(num_seqs):
            for f in range(num_feats):
                raw_val = compfeatvals_all[f][s]
                if raw_val == 0 or stdvals_species[f] < min_std_threshold:
                    subzveccomp[s, f] = np.nan
                else:
                    subzveccomp[s, f] = (raw_val - meanvals_species[f]) / stdvals_species[f]

        raw_cols = [f"nardini_{feat}_raw" for feat in mycompfeats_all]
        zscore_cols = [f"nardini_{feat}_zscore" for feat in mycompfeats_all]
        raw_df.columns = raw_cols
        zscore_df = pd.DataFrame(data=subzveccomp, columns=zscore_cols)
        nardini_df = pd.concat([raw_df, zscore_df], axis=1)
        print(f"Calculated {num_feats} raw features and {num_feats} z-scores.\n")
    else:
        raw_cols = [f"nardini_{feat}_raw" for feat in mycompfeats_all]
        raw_df.columns = raw_cols
        nardini_df = raw_df
        print(f"Calculated {num_feats} raw features (no z-scores).\n")

    # ========================================================================
    # CALCULATE CIDER AND SPARROW PARAMETERS
    # ========================================================================
    print("Calculating CIDER and SPARROW parameters...")

    cider_data = []
    sparrow_data = []

    for idx, currseq in enumerate(subseqs):
        if len(currseq) < 1 or any(x in currseq for x in "XUZJBO"):
            cider_data.append({'cider_kappa': np.nan, 'cider_omega': np.nan, 'cider_delta': np.nan,
                                'cider_uversky_hydropathy': np.nan, 'cider_fraction_neutral': np.nan})
            sparrow_data.append({'sparrow_SCD': np.nan, 'sparrow_SHD': np.nan, 'sparrow_complexity': np.nan,
                                  'sparrow_fraction_proline': np.nan, 'sparrow_scaled_rg': np.nan,
                                  'sparrow_scaled_re': np.nan, 'sparrow_prefactor': np.nan,
                                  'sparrow_scaling_exponent': np.nan, 'sparrow_asphericity': np.nan})
            continue

        try:
            SeqOb = CIDER_SP(currseq)
            cider_data.append({
                'cider_kappa': SeqOb.get_kappa(),
                'cider_omega': SeqOb.get_Omega(),
                'cider_delta': SeqOb.get_delta(),
                'cider_uversky_hydropathy': SeqOb.get_uversky_hydropathy(),
                'cider_fraction_neutral': SeqOb.get_countNeut() / SeqOb.get_length(),
                'cider_length': SeqOb.get_length()
            })
        except Exception as e:
            print(f"  Warning: CIDER failed for sequence {idx+1}: {e}")
            cider_data.append({'cider_kappa': np.nan, 'cider_omega': np.nan, 'cider_delta': np.nan,
                                'cider_uversky_hydropathy': np.nan, 'cider_fraction_neutral': np.nan,
                                'cider_length': np.nan})

        try:
            prot = Protein(currseq)
            sparrow_dict = {
                'sparrow_SCD': prot.SCD,
                'sparrow_SHD': prot.SHD,
                'sparrow_complexity': prot.complexity,
                'sparrow_fraction_proline': prot.fraction_proline
            }
            try:
                sparrow_dict['sparrow_scaled_rg'] = prot.predictor.radius_of_gyration(use_scaled=True)
            except Exception:
                sparrow_dict['sparrow_scaled_rg'] = np.nan
            try:
                sparrow_dict['sparrow_scaled_re'] = prot.predictor.end_to_end_distance(use_scaled=True)
            except Exception:
                sparrow_dict['sparrow_scaled_re'] = np.nan
            try:
                sparrow_dict['sparrow_prefactor'] = prot.predictor.prefactor()
            except Exception:
                sparrow_dict['sparrow_prefactor'] = np.nan
            try:
                sparrow_dict['sparrow_scaling_exponent'] = prot.predictor.scaling_exponent()
            except Exception:
                sparrow_dict['sparrow_scaling_exponent'] = np.nan
            try:
                sparrow_dict['sparrow_asphericity'] = prot.predictor.asphericity()
            except Exception:
                sparrow_dict['sparrow_asphericity'] = np.nan

            sparrow_data.append(sparrow_dict)
        except Exception as e:
            print(f"  Warning: SPARROW failed for sequence {idx+1}: {e}")
            sparrow_data.append({'sparrow_SCD': np.nan, 'sparrow_SHD': np.nan, 'sparrow_complexity': np.nan,
                                  'sparrow_fraction_proline': np.nan, 'sparrow_scaled_rg': np.nan,
                                  'sparrow_scaled_re': np.nan, 'sparrow_prefactor': np.nan,
                                  'sparrow_scaling_exponent': np.nan, 'sparrow_asphericity': np.nan})

    cider_df = pd.DataFrame(cider_data)
    sparrow_df = pd.DataFrame(sparrow_data)

    print("CIDER and SPARROW calculation complete.\n")

    # ========================================================================
    # COMBINE ALL FEATURES AND SAVE
    # ========================================================================
    name_df = pd.DataFrame({'Name': subnames})
    seq_df = pd.DataFrame({'Sequence': subseqs})
    final_df = pd.concat([name_df, seq_df, nardini_df, cider_df, sparrow_df], axis=1)

    print("="*70)
    print("CHUNK SUMMARY")
    print(f"Sequences in this chunk: {len(final_df)}")
    print(f"Total columns: {len(final_df.columns)}")
    print("="*70 + "\n")

    final_df.to_csv(args.output_csv, index=False)
    print(f"Saved chunk output to: {args.output_csv}")

if __name__ == "__main__":
    main()
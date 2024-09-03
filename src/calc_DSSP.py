# -*- coding: utf-8 -*-
#!/usr/bin/python
"""
@author: jmrodriguezc
"""

# import global modules
import sys
import os
import argparse
import logging
import re
import glob
import pandas as pd
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.SeqUtils import seq1
import concurrent.futures
from itertools import repeat



###################
# Parse arguments #
###################

parser = argparse.ArgumentParser(
    description='Calculate the RSA from DSSP',
    epilog='''Examples:
        
    python  calc_RSA.py
        -w 10
        -i databases/AlphaFold/202403_F1-model_v4/homo_sapiens
        -o databases/DSSP/202403_F1-model_v4/homo_sapiens
    ''',
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-w',   type=int, default=4, help='Number of threads/n_workers (default: %(default)s)')
parser.add_argument('-i',   required=True, help='Input dir with the PDB files')
parser.add_argument('-o',   required=True, help='Output dir for the DSSP results')
args = parser.parse_args()

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')


#############
# Constants #
#############
AF_PDB_PATTHER = 'AF-(.*)-F1-model_v4'

DSSP_HEADER = [
'DSSP index',
'Amino acid',
'Secondary structure',
'Relative ASA',
'Phi',
'Psi',
'NH->O_1_relidx',
'NH->O_1_energy',
'O->NH_1_relidx',
'O->NH_1_energy',
'NH->O_2_relidx',
'NH->O_2_energy',
'O->NH_2_relidx',
'O->NH_2_energy']


###################
# Local functions #
###################
# Calculate the 
# rsa_window (default 25) - RSA values are smoothed over a window centered on the residue to predict
# rsa_threshold (default 0.581) - Binding predictions are overweighted when disorder prediction is above this threshold
# 
# https://github.com/BioComputingUP/AlphaFold-disorder/
# 
# Piovesan D, Monzon AM, Tosatto SCE.
# Intrinsic protein disorder and conditional folding in AlphaFoldDB. Protein Sci. 2022 Nov;31(11):e4466.
# PMID: 36210722 PMCID: PMC9601767.

def moving_average(x, w):
    # https://stackoverflow.com/questions/13728392/moving-average-or-running-mean
    return np.convolve(x, np.ones(w), 'valid') / w

def make_prediction(df, window_rsa=[25], thresholds_rsa=[0.581]):
    for w in window_rsa:
        # Smooth disorder score (moving average)
        column_rsa_window = 'disorder-{}'.format(w)
        half_w = int((w - 1) / 2)
        df[column_rsa_window] = moving_average(np.pad(df['Relative ASA'], (half_w, half_w), 'reflect'), half_w * 2 + 1)

        # Transofrm scores above RSA threshold
        for th_rsa in thresholds_rsa:
            column_rsa_binding = 'binding-{}-{}'.format(w, th_rsa)
            df[column_rsa_binding] = df[column_rsa_window].copy()
            df.loc[df[column_rsa_window] > th_rsa, column_rsa_binding] = df.loc[
                                            df[column_rsa_window] > th_rsa, 'pLDDT'] * (1 - th_rsa) + th_rsa
    return df

# Calculate the DSSP from protein id
def calc_dssp(pdb, odir):
    # init variables
    dssp_report = []
    df = pd.DataFrame()
    # get the file name (pdb name)
    pdb_code = os.path.splitext(os.path.basename(pdb))[0]
    # get the protein id
    match = re.search(AF_PDB_PATTHER, pdb_code)
    q_id = match.group(1) if match else pdb_code
    # calculate the DSSP
    try:
        # parse the PDB
        p = PDBParser()
        structure = p.get_structure(q_id, pdb)
        # model the PDB
        model = structure[0]
        dssp = DSSP(model, pdb, dssp='mkdssp')
        # get the keys from DSSP
        dssp_report = [dssp[k] for k in list(dssp.keys())]
    except Exception as e:
        logging.error(f"ERROR: calculating the DSSP for {q_id}: {e}")
    
    # create dataframe for DSSP
    try:
        if len(dssp_report) > 0:
            # create and print df
            df = pd.DataFrame(list(dssp_report), columns=DSSP_HEADER)
    except Exception as e:
        logging.error(f"ERROR: creating the dataframe for {q_id}: {e}")

    if not df.empty:
        # parse b-factor (pLDDT) and DSSP
        try:
            df2 = []
            for i, residue in enumerate(model.get_residues()):
                lddt = residue['CA'].get_bfactor()
                df2.append((i + 1, seq1(residue.get_resname()), lddt))
            df2 = pd.DataFrame(df2, columns=['pos', 'aa', 'pLDDT'])
            # merge the DSSP and pLDDT scores
            df = pd.merge(df, df2, left_on=['DSSP index','Amino acid'], right_on=['pos','aa'])
            # remove columns
            df.drop(['pos','aa'], axis=1, inplace=True)
        except Exception as e:
            logging.error(f"ERROR: parsing b-factor (pLDDT) and RSA for {q_id}: {e}")


        # make prediction:
        # rsa_window - RSA values are smoothed over a window centered on the residue to predict
        # rsa_threshold - Binding predictions are overweighted when disorder prediction is above this threshold
        try:
            df = make_prediction(df.copy())
        except Exception as e:
            logging.error(f"ERROR: making smoothed over a window for {q_id}: {e}")

    
        # printing DSSP+AlphaFold report
        try:
            ofile = f"{odir}/DSSP_{pdb_code}.tsv"
            df.to_csv(ofile, sep="\t", index=False)
        except Exception as e:
            logging.error(f"ERROR: printing DSSP+AlphaFold for {q_id}: {e}")
       
        # add the protein id into dataframe
        try:
            df.insert(0, 'Protein_ID', q_id)
        except Exception as e:
            logging.error(f"ERROR: adding the protein id into dataframe for {q_id}: {e}")

    return df
    

#################
# Main function #
#################
def main(args):
    '''
    Main function
    '''
    
    logging.info("getting the input parameters...")
    n_workers = args.w
    idir = args.i
    odir = args.o
    ofile = f"{odir}/DSSP.tsv"


    logging.info("preparing the workspace...")
    os.makedirs(odir, exist_ok=True)
    
    
    logging.info("getting the PDB files...")
    ipdbs = glob.glob(f"{idir}/*.pdb")
    
    
    logging.info("calculating the DSP for all given PDBs...")
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as executor:        
        dssp_results = executor.map( calc_dssp, ipdbs, repeat(odir))
    dssp_results = pd.concat(dssp_results)
    # begin:
    # dssp_results = calc_dssp(ipdbs[0], odir)
    # end


    logging.info("printing the global DSSP report...")
    dssp_results.to_csv(ofile, sep="\t", index=False)
        
        
    logging.info("the program has finished sucessfully")


if __name__ == "__main__":
    # start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    main(args)
    logging.info('end script')


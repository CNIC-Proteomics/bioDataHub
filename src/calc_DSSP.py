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
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
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
# Calculate the DSSP from protein id
def calc_dssp(pdb, odir):
    # init variables
    dssp_report = []
    df = pd.DataFrame()
    # get the file name (pdb name)
    pdb_name = os.path.splitext(os.path.basename(pdb))[0]
    # get the protein id
    match = re.search(AF_PDB_PATTHER, pdb_name)
    q_id = match.group(1) if match else pdb_name
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
    
    # create dataframe and print DSSP
    try:
        if len(dssp_report) > 0:
            # create and print df
            df = pd.DataFrame(list(dssp_report), columns=DSSP_HEADER)
            ofile = f"{odir}/DSSP_{pdb_name}.tsv"
            df.to_csv(ofile, sep="\t", index=False)
            # add the protein id into dataframe
            df.insert(0, 'Protein_ID', q_id)
    except Exception as e:
        logging.error(f"ERROR: creating the dataframe for {q_id}: {e}")

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
    # begin:
    # dssp_results = calc_dssp(ipdbs[0], odir)
    # end
    dssp_results = pd.concat(dssp_results)


    logging.info("printing the global DSSP report...")
    dssp_results.to_csv(ofile, sep="\t", index=False)
        
        
    logging.info("the program has finished sucessfully")


if __name__ == "__main__":
    # start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    main(args)
    logging.info('end script')


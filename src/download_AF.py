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
import pandas as pd
import numpy as np
import urllib.request
import json
import glob
import tarfile
import gzip
import concurrent.futures
from itertools import repeat



###################
# Parse arguments #
###################

parser = argparse.ArgumentParser(
    description='Download the AlphaFold proteomes',
    epilog='''Examples:
        
    python download_AF.py
        -w 10
        -i filtered_LIMMA_NM_pgmqfall_table.APPRIS.tsv
        -o DSSP
    ''',
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-w',   type=int, default=4, help='Number of threads/n_workers (default: %(default)s)')
parser.add_argument('-s',   required=True, help='Species name')
parser.add_argument('-v',   required=True, help='AlphaFold version')
parser.add_argument('-f',   required=True, help='FASTA file to get the proteins')
parser.add_argument('-o',   required=True, help='Output folder')
args = parser.parse_args()

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')


#############
# Constants #
#############

# AlphaFold constants
AF_EBI_URL = 'https://alphafold.ebi.ac.uk/files/'
AF_PDB_NAME = 'AF-__PROTEINID__-F1-model_--VERSION--.pdb'

# Local folder
LOCAL_DIR = os.path.dirname(os.path.abspath(__file__))

# Config file for the species
with open(os.path.join(LOCAL_DIR, 'config.json')) as f:
   SPECIES_CFG = json.load(f)

###################
# Local functions #
###################

# Parse FASTA file
def parse_fasta(fasta_file):
    sequences = {}
    current_sequence = None
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # New sequence header
                current_sequence = line[1:].split()[0]
                sequences[current_sequence] = ''
            elif current_sequence is not None:
                # Append sequence data
                sequences[current_sequence] += line
    return sequences


# Extract the UniProt accesions
def extract_uniprot_accessions(fasta_file):
    sequences = parse_fasta(fasta_file)
    accessions = [header.split('|')[1] for header in sequences.keys()]
    return accessions


# Uncompress tar file
def extract_tar_file(ifile, odir):
    try:
        with tarfile.open(ifile, 'r') as tar:
            # Extract all contents to the destination directory
            tar.extractall(odir)
    except Exception as exc:
        logging.error(f"ERROR: failed uncompressing the file {ifile}: {exc}")


# Uncompress "gz" file within folder
def uncompress_gz_files(ifile, odir):
    try:
        # Create the output filename by removing the ".gz" extension
        ofile = os.path.splitext(ifile)[0]
        # Read the gzipped file and write the uncompressed content
        with gzip.open(os.path.join(odir, ifile), 'rb') as f_in:
            with open(os.path.join(odir, ofile), 'wb') as f_out:
                f_out.write(f_in.read())
    except Exception as exc:
        logging.warning(f"failed uncompressing the file {ofile}: {exc}")
    # return the ibase file name (without .gz extension)
    return os.path.basename( ofile )


# Download AF proteome
def dowloand_af(url, odir):
    
    # extract filename from the URL
    ofile = url.split('/')[-1]
    ofile = ofile.split('?')[0] # if there's a query string, remove it
    ofile = os.path.join(odir,ofile)
    
    # download the file
    logging.info(f"wget {url} > {ofile}")
    try:
        urllib.request.urlretrieve(url, ofile)
    except Exception as exc:
        logging.error(f"ERROR: failed to dowload {url}: {exc}")
    
    # Uncompress the tar file
    logging.info("uncompressing AlphaFold tar file...")
    extract_tar_file(ofile, odir)
    
    # return the compressed files
    ofiles = glob.glob(f"{odir}/*.pdb.gz")
    return ofiles

    

# Try to download more PDB files if there is not in the AF proteome
def dowloand_more_af(qid, pdb_files, idir):
    
    # get the name of PDF in AF
    pdb_name = AF_PDB_NAME.replace('__PROTEINID__', qid)
    
    # check if exists the loal PDB; otherwise, download from AF-EBI
    if not pdb_name in pdb_files:
        try:
            # download the AF PDB from EBI
            url = AF_EBI_URL+pdb_name
            ofile = f"{idir}/{pdb_name}"
            logging.info(f"wget {url} > {ofile}")
            urllib.request.urlretrieve(url, ofile)
        except Exception as e:
            logging.warning(f"downloading from EBI web the AF for {qid}: {e}")

    

#################
# Main function #
#################
def main(args):
    '''
    Main function
    '''
    
    logging.info("getting the input parameters...")
    n_workers = args.w
    species = args.s
    version = args.v
    fastafile = args.f
    odir = args.o
    # species = "human"
    # version = "v4"
    # fastafile = r"S:\U_Proteomica/UNIDAD/Databases/UniProt/202311/human/sequences/human_202311_pro-sw-tr.fasta"
    # odir = r"S:\U_Proteomica/UNIDAD/Databases/AlphaFold/202403.v4/test"


    logging.info("preparing the workspace...")
    os.makedirs(odir, exist_ok=True)
    


    logging.info("dowloading the AlphaFold proteome...")
    if 'alphafold' in SPECIES_CFG[species]:
        AF_PROTEOME_URL = SPECIES_CFG[species]['alphafold']
        AF_PROTEOME_URL = AF_PROTEOME_URL.replace('--VERSION--',version)
        global AF_PDB_NAME
        AF_PDB_NAME = AF_PDB_NAME.replace('--VERSION--',version)
    else:
        logging.error("ERROR: the species has not AlphaFold proteome repository")
        sys.exit()
    pdbgz_files = dowloand_af(AF_PROTEOME_URL, odir)



    logging.info("uncompressing the PDB files from AlphaFold...")
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as executor:
        pdb_files = executor.map( uncompress_gz_files, pdbgz_files, repeat(odir) )
    pdb_files = list(pdb_files)
    # begin:
    # for Spyder
    # pdb_files = uncompress_gz_files(pdbgz_files[0], odir)
    # end

    
    logging.info("extracting the UniProt accessions from FASTA file...")
    qs = extract_uniprot_accessions(fastafile)
    qs = np.unique(qs)


   
    logging.info("dowloading more AlphaFold predictions from the given proteins...")
    with concurrent.futures.ProcessPoolExecutor(max_workers=n_workers) as executor:
        executor.map( dowloand_more_af, qs, repeat(pdb_files), repeat(odir) )
    # begin:
    # for Spyder
    # dowloand_more_af(qs[0], pdb_files, odir)
    # end


    
    logging.info("the program has finished sucessfully")


if __name__ == "__main__":
    # start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    main(args)
    logging.info('end script')


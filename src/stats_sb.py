#!/usr/bin/python
import sys
import argparse
import logging
import pandas as pd
import numpy as np

__author__ = 'jmrodriguezc'
__credits__ = ["Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "2.4"
__maintainer__ = "Jose Rodriguez"
__email__ = "jmrodriguezc@cnic.es"
__status__ = "Development"


def main(args):
    ''' Main function'''

    logging.info("openning category file...")
    cat = pd.read_csv(args.infile, sep="\t", low_memory=False)

    logging.info("getting the number of reviewed proteins and unreviewed proteins with isoforms!...")
    # get the number of reviewed proteins and unreviewed proteins (with isoforms!!!)
    stats = cat[['Protein','prot_UniProt_Class']]    
    reviewed_isof = stats.groupby(['prot_UniProt_Class'])['Protein'].nunique().reset_index()
    reviewed_isof.columns = range(reviewed_isof.columns.size) # reset columns
    # # get the number of reviewed proteins and unreviewed proteins (with isoforms!!!)
    # stats = stats[~stats['Protein'].str.contains('-')]
    # reviewed_prot = stats.groupby(['prot_UniProt_Class'])['Protein'].nunique().reset_index()
    # reviewed_prot.columns = range(reviewed_prot.columns.size) # reset columns

    logging.info("getting the number of isoforms per category...")
    # get the number of isoforms per category
    # add the column name to all values, except is null
    # create a dict with the number of proteins per category
    c_cols = [ c for c in cat.columns if c.startswith('cat_') or c == 'APPRIS Annotation' or c == 'norm_trifid_score' or c == 'corsair_score']
    c_prot = cat[['Protein']+c_cols]
    out = dict()
    for c in c_prot.columns:
        if c in c_cols:
            c_prot.loc[c_prot[c].notnull(), c] = c    
            out.update( c_prot.groupby(c)['Protein'].nunique().to_dict() )
    category_prot = pd.DataFrame().from_dict(out.items())

    logging.info("concating the stats reports...")
    out = pd.concat([reviewed_isof,category_prot])
    out.columns = ['Label','Isoforms']

    logging.info('printing stats file...')
    out.to_csv(args.outfile, index=False, sep="\t")



if __name__ == "__main__":
    # parse arguments
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass
    parser = argparse.ArgumentParser(
        description='Get the statistics from category file',
        epilog='''
        ''',
        formatter_class=CustomFormatter )
    parser.add_argument('-i',  '--infile', required=True, help='Input file')
    parser.add_argument('-o',  '--outfile', required=True, help='Output file')
    parser.add_argument('-v', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()

    # logging debug level. By default, info level
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(level=logging.INFO,
                            format='%(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')

    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    main(args)
    logging.info('end script')

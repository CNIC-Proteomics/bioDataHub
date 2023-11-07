# -*- coding: utf-8 -*-
#!/usr/bin/python

# Module metadata variables
__author__ = ["Jose Rodriguez"]
__credits__ = ["Jose Rodriguez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.0.1"
__maintainer__ = "Jose Rodriguez"
__email__ = "jmrodriguezc@cnic.es"
__status__ = "Development"

# import global modules
import sys
import argparse
import logging
import pandas as pd
import numpy as np

###################
# Parse arguments #
###################

parser = argparse.ArgumentParser(
    description='Extract categories from orthologous genes',
    epilog='''Examples:
        
    python  categorize_orthologous.py
      -im  data/rabbit/mart_export.tsv
      -ic1 data/human_202306.uniprot.tsv
      -ic2 data/rabbit/rabbit_202306.uniprot.tsv
      -o   data/rabbit/rabbit_202306_with_human_orthologs.uniprot.tsv
    ''',
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-im',   required=True, help='BioMart results showing the relationship between the human gene ortholog and genes in other species')
parser.add_argument('-ic1',  required=True, help='Human categories from UniProtKB')
parser.add_argument('-ic2',  required=True, help='Categories from other species in UniProtKB')
parser.add_argument('-o',    required=True, help='Results with the categories from human orthologous')

args = parser.parse_args()

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')


#################
# Main function #
#################
def main(args):
    '''
    Main function
    '''    
    logging.info("getting the input parameters...")
    ifile_m = args.im
    ifile_c1 = args.ic1
    ifile_c2 = args.ic2
    ofile = args.o


    logging.info("reading input files...")
    data_m = pd.read_csv(ifile_m, sep="\t", low_memory=False)
    data_c1 = pd.read_csv(ifile_c1, sep="\t", low_memory=False)
    data_c2 = pd.read_csv(ifile_c2, sep="\t", low_memory=False)
    cats = [c for c in data_c1.columns if c.startswith('cat_')]

    logging.info("filtering the orthologs by high confidence score...")
    data_m = data_m[data_m['Human orthology confidence [0 low, 1 high]'] == 1]
    # data_m = data_m[['Gene stable ID','Transcript stable ID','Protein stable ID','Human protein or transcript stable ID']]
    data_m = data_m[['Protein stable ID','Human protein or transcript stable ID']].drop_duplicates()
    data_m = data_m.dropna() # drop rows woth any nan

    logging.info("obtaining the list of human orthologous proteins...")
    orthologs_human = np.unique(data_m['Human protein or transcript stable ID'].tolist())
    

    # logging.info("obtaining the list of human orthologous proteins...")
    # orthologs_species = np.unique(data_m['Protein stable ID'].tolist())

    
    logging.info("filtering the human categories based on the ensembl protein/transcript ids...")
    data_human = data_c1[data_c1['xref_Ensembl_protId'].isin(orthologs_human) | data_c1['xref_Ensembl_transcId'].isin(orthologs_human)]
    # data_human = data_c1.merge(data_m2, left_on=['xref_Ensembl_protId'], right_on=['Human protein or transcript stable ID'])
    data_human = data_human[['xref_Ensembl_protId']+cats]
    data_human = data_human.drop_duplicates()
    


    logging.info("filtering the species report based on protein ID...")
    # data_species = data_c2[data_c2['xref_Ensembl_protId'].isin(orthologs_species)]
    # data_species2 = data_c2.merge(data_m2, left_on=['xref_Ensembl_protId'], right_on=['Protein stable ID'])
    data_species = data_c2.merge(data_m, how='left', left_on=['xref_Ensembl_protId'], right_on=['Protein stable ID']).drop_duplicates()
    
    
    logging.info("Adding the human categories to the other species if they do not exist...")
    # merge both datasets
    data_species = data_species.merge(data_human, how='left', left_on=['Human protein or transcript stable ID'], right_on=['xref_Ensembl_protId']).drop_duplicates()
    # fill the nan categories with the human categories
    for c in cats:
        data_species[c] = data_species[f"{c}_x"].fillna(data_species[f"{c}_y"])
    # remove the temporal categories (suffixes 'x' and 'y')
    cc = [c for c in data_species.columns if c.startswith('cat_') and (c.endswith('_x') or c.endswith('_y')) or c == 'xref_Ensembl_protId_y' ]
    data_species = data_species.drop(columns=cc)
    # remove duplicates
    data_species = data_species.drop_duplicates()


    logging.info("printing the output file...")
    data_species.to_csv(ofile, sep="\t", index=False)


if __name__ == "__main__":
    # start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    main(args)
    logging.info('end script')


# 1. 
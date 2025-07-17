#!/usr/bin/python
__author__ = 'jmrodriguezc'
__credits__ = ["Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.0.2"
__maintainer__ = "Jose Rodriguez"
__email__ = "jmrodriguezc@cnic.es"
__status__ = "Development"

import sys
import os
import argparse
import logging
import pandas as pd
import glob


####################
# Global variables #
####################
GTF_COLUMNS = [
    'seqname',
    'source',
    'feature',
    'start',
    'end',
    'score',
    'strand',
    'frame',
    'attribute'
]

OUTPUT_COLUMNS = [
    'seqname',
    'source',
    'feature',
    'pep_start',
    'pep_end',
    'score',
    'strand',
    'frame',
    'ensembl_gene_id',
    'ensembl_transc_id',
    'note'
]

REF_COLUMNS = ['seqname','pep_start','pep_end']

###################
# Local functions #
###################
def convert_to_protein_gtf(file):

    logging.info(f"{os.path.basename(file)}")
    df = pd.read_csv(file, sep="\t", dtype=str, header=None, names=GTF_COLUMNS, low_memory=False)
    # df_bak = pd.read_csv(file, sep="\t", dtype=str, header=None, names=GTF_COLUMNS, low_memory=False)
    # df = df_bak.head(1000)
    
    # get the name of last column (attributes)
    col_attr = df.columns[-1]

    # create a new df with the note values in separate columns, if they exist
    df_notes = pd.DataFrame()
    df_notes['ensembl_gene_id'] = df[col_attr].str.extract(r'gene_id \"([^\"]*)\"', expand=False)
    df_notes['ensembl_transc_id']= df[col_attr].str.extract(r'transcript_id \"([^\"]*)\"', expand=False)
    flag_notes = False
    if df[col_attr].str.contains('note').any(): # get the position for appris methods
        df_notes['note']= df[col_attr].str.extract(r'note \"([^\"]*)\"', expand=False)
        # check if pep_start and pep_end exist
        if df[col_attr].str.contains('pep_start').any() and df[col_attr].str.contains('pep_end').any():
            # extract info
            df_notes['pep_start']= df[col_attr].str.extract(r'pep_start:(\d+)', expand=False)
            df_notes['pep_end']= df[col_attr].str.extract(r'pep_end:(\d+)', expand=False)
            # remove duplicated information in note column
            df_notes['note'] = df_notes['note'].str.replace(',*pep_start:\d*,*','',regex=True).str.replace(',*pep_end:\d*,*','',regex=True)
            flag_notes = True
        elif df[col_attr].str.contains('pep_position').any():
            # extract info
            df_notes['pep_start']= df[col_attr].str.extract(r'pep_position:(\d+)', expand=False)
            df_notes['pep_end']= df_notes['pep_start']
            # remove duplicated information in note column
            df_notes['note'] = df_notes['note'].str.replace(',*pep_position:\d*,*','',regex=True)
            flag_notes = True
    if df[col_attr].str.contains('reliability').any(): # get aa length for the appris annotation. add appris label in the reliability attribute
        # add appris label
        df_notes['note']= df[col_attr].str.extract(r'reliability \"([^\"]*)\"', expand=False)
        # extract info
        df_notes['pep_start']= '1'
        df_notes['pep_end']= df[col_attr].str.extract(r'length_aa \"([^\"]*)\"', expand=False)
        flag_notes = True

    # if the "df notes" is not empty, create a new df with GTF format for proteins
    df_gtf = pd.DataFrame()
    if flag_notes:
        df_gtf = pd.concat([df,df_notes], axis=1)
        # select the specific columns
        df_gtf = df_gtf[OUTPUT_COLUMNS]
        # update the values of seqname with ensembl_transc_id
        df_gtf['seqname'] = df_gtf['ensembl_transc_id']
        # drop the rows that seqname,pep_start, or pep_end are empty
        df_gtf = df_gtf.dropna(subset = REF_COLUMNS)
           
    return df_gtf

#################
# Main function #
#################
def main(args):
    '''
    Main function
    '''
    
    logging.info("processing the identification files...")
    appris_files = glob.glob(args.input_appris)

    
    logging.info("converting to protein report in GTF format...")
    appris_gtfs = [ convert_to_protein_gtf(f) for f in appris_files ]

    
    logging.info("merging the whole protein reports based on APPRIS methods...")
    appris_gtf = pd.concat(appris_gtfs)
    # sort by REF_COLUMNS
    appris_gtf = appris_gtf.sort_values(by=REF_COLUMNS)
    
    
    logging.info("preparing the uniprot-ensembl relationship...")
    uniprot_df = pd.read_csv(args.input_uniprot, sep="\t", dtype=str, na_values=['NA'], low_memory=False)
    # get the relationship UniProt-Ensembl (transcript id)
    uniprot_df = uniprot_df[['Gene','Protein','xref_Ensembl_transcId']]
    uniprot_df = uniprot_df.dropna()
    # rename tables
    uniprot_df.rename(columns={'Protein': 'uniprot_id', 'Gene': 'gene_name'}, inplace=True)


    logging.info("merging the relationship UniProt-Ensembl(transc_id)...")
    appris_gtf = appris_gtf.merge(uniprot_df, how='left', left_on=['ensembl_transc_id'], right_on=['xref_Ensembl_transcId']).drop_duplicates()
    # remove obsolete cols
    appris_gtf = appris_gtf.drop(columns='xref_Ensembl_transcId')
    
    
    logging.info(f"preparing the {args.type_out} output...")
    if args.type_out == 'uniprot':
        appris_gtf['seqname'] = appris_gtf['uniprot_id']
    # drop the rows that seqname,pep_start, or pep_end are empty
    appris_gtf = appris_gtf.dropna(subset = REF_COLUMNS)
    # add the slash in the header (in the first colum name)
    appris_gtf.rename(columns={'seqname': '#seqname'}, inplace=True)
    
    
    logging.info("printing output...")
    appris_gtf.to_csv(args.outfile, sep="\t", index=False)
    
    
    



if __name__ == "__main__":
    # parse arguments
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass
    parser = argparse.ArgumentParser(
        description='Convert the method annotations in GTF format to another GTF that references the protein region',
        epilog='''
Examples:
    python  src/convert_appris.py  -ia test/appris/*.gtf  -iu test/human_202306.uniprot.tsv -o test/appris/human_202306.appris.gtf
        ''',
        formatter_class=CustomFormatter )
    parser.add_argument('-ia',  '--input-appris', required=True, help='APPRIS annotation file(s) in GTF format')
    parser.add_argument('-iu',  '--input-uniprot', required=True, help='UniProt report with the cross-reference identifiers')
    parser.add_argument('-t',   '--type-out', choices=['uniprot','ensembl'], default='uniprot', help='Determine the type of database to reference in the output file')
    parser.add_argument('-o',   '--outfile', required=True, help='Output file')
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

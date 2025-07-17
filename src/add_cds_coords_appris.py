#!/usr/bin/python
__author__ = 'jmrodriguezc'
__credits__ = ["Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.0.2"
__maintainer__ = "Jose Rodriguez"
__email__ = "jmrodriguezc@cnic.es"
__status__ = "Development"

import sys
import argparse
import logging
import pandas as pd


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
    'attributes'
]

PROT_COLUMNS = [
    'seqname',
    'source',
    'feature',
    'pep_start',
    'pep_end',
    'score',
    'strand',
    'frame',
    'attributes'
]

###################
# Local functions #
###################


#################
# Main function #
#################
def main(args):
    '''
    Main function
    '''
    
    logging.info("getting the input parameters...")
    ia_ifile = args.input_appris
    ip_ifile = args.input_cds
    # ofile    = args.output
    ic_ofile = args.output_cds
    # ia_ifile = r"S:\U_Proteomica\UNIDAD\Softwares\jmrodriguezc\SANPRO\tmp\202501\zebrafish\appris_method.appris.gtf"
    # ip_ifile = r"S:\U_Proteomica\UNIDAD\Softwares\jmrodriguezc\SANPRO\tmp\202501\zebrafish\appris_data.pannot.gtf"
    # ofile    = r"S:\U_Proteomica\UNIDAD\Softwares\jmrodriguezc\SANPRO\data\202501\zebrafish\zebrafish_202501.appris.cds.tsv"
    # ic_ofile = r"S:\U_Proteomica\UNIDAD\Softwares\jmrodriguezc\SANPRO\data\202501\zebrafish\zebrafish_202501.appris.cds.gtf"



    logging.info("reading APPRIS annotation table...")
    appris_gtf = pd.read_csv(ia_ifile, sep="\t", dtype=str, header=None, low_memory=False, names=GTF_COLUMNS)
    # filter the APPRIS annot
    appris_gtf = appris_gtf[appris_gtf['feature'] == 'principal_isoform']



    logging.info("extracting APPRIS annotation...")
    # extract transcript_id from the attributes column
    attribute_regex = r'.*transcript_id\s"([^"]+)"'
    appris_gtf[['transcript_id']] = appris_gtf['attributes'].str.extract(attribute_regex)
    # extract gene_id from the attributes column
    attribute_regex = r'.*gene_id\s"([^"]+)"'
    appris_gtf[['gene_id']] = appris_gtf['attributes'].str.extract(attribute_regex)
    # extract reliability from the attributes column
    attribute_regex = r'.*reliability\s"([^"]+)"'
    appris_gtf[['reliability']] = appris_gtf['attributes'].str.extract(attribute_regex)
    
    
    
    logging.info("reading Protein table with the CDS coordinates...")
    pannot_df = pd.read_csv(ip_ifile, sep="\t", header=None, low_memory=False, names=PROT_COLUMNS)
    pannot_df = pannot_df.astype(str)
    # remove bad strand coming from protein annotation
    pannot_df = pannot_df.drop(columns=['strand'], axis=1)




    logging.info("parsing the attributes column of protein table...")
    # split the attributes column by ";" into multiple columns
    pannot_attrs = pannot_df['attributes'].str.split(';', expand=True)
    # remove labels (like ID=, Parent=, Gene=, Note=) and create columns with their names
    for col in pannot_attrs.columns:
        if pannot_attrs[col].notnull().all():
            # Extract the key (e.g., ID, Parent, Gene) and value
            pannot_attrs[col] = pannot_attrs[col].str.strip()
            if "=" in pannot_attrs[col].iloc[0]:
                key = pannot_attrs[col].iloc[0].split('=')[0].strip()
                pannot_attrs[col] = pannot_attrs[col].str.split('=', expand=True)[1]
                pannot_attrs = pannot_attrs.rename(columns={col: key})
    # rename columns
    pannot_attrs = pannot_attrs.rename(columns={'ID': 'exon_id', 'Parent': 'transcript_id', 'Gene': 'gene_id', 'Note': 'cds_coords'})
    # remove the substring
    pannot_attrs['cds_coords'] = pannot_attrs['cds_coords'].str.replace('cds_coord>', '')
    # get the cds coordinates in separated columns
    pannot_attrs[['coords','strand']] = pannot_attrs['cds_coords'].str.split(':', expand=True)
    pannot_attrs[['cds_start','cds_end']] = pannot_attrs['coords'].str.split('-', expand=True)
    pannot_attrs = pannot_attrs.drop(columns=['cds_coords','coords'], axis=1)
    # combine the split attributes back with the original DataFrame
    pannot_df = pd.concat([pannot_df.drop(columns=['attributes']), pannot_attrs], axis=1)




    logging.info("adding the chr value into CDS coordinate of protein table...")
    # get 'seqname' (chr) and transcript_id without duplicates
    appris_chr_df = appris_gtf[['transcript_id','seqname','reliability']].drop_duplicates().rename(columns={'seqname':'chr', 'reliability': 'appris_label'})
    # merge the two DataFrames on gene_id. We'll add the seqname from appris_gtf to pannot_df
    pannot_df = pannot_df.merge(appris_chr_df, on='transcript_id', how='left')
    # rename feature attribute with the appris annotation
    pannot_df['feature'] = pannot_df['appris_label']
    pannot_df['source'] = "APPRIS"
    # get columns
    pannot_df = pannot_df[[
        'source','feature','pep_start','pep_end',
        'chr','cds_start','cds_end','strand','frame',
        'exon_id', 'transcript_id', 'gene_id',
        'score'
    ]]
    
    
    

    
    logging.info("create ouputs in GTF format...")    
    # create the report with the CDS coordinates
    # add the attributes column in GTF format with the rest of attributes
    cannot_gtf = pannot_df[[
        'chr',
        'source',
        'feature',
        'cds_start',
        'cds_end',
        'score',
        'strand',
        'frame'
    ]].copy()
    # get attributes
    cannot_gtf['attributes'] = ""
    cc = ['exon_id','transcript_id','gene_id','chr','cds_start','cds_end','strand','frame','pep_start','pep_end','feature']
    for c in cc:
        if c == 'feature':
            cannot_gtf['attributes'] += 'appris_label=' + pannot_df[c] + ';'
        else:    
            cannot_gtf['attributes'] += c + '=' + pannot_df[c] + ';'
    # overwrite values
    cannot_gtf['feature'] = "CDS"
    cannot_gtf['source'] = "GENCODE"



    logging.info("printing output...")
    # pannot_df.to_csv(ofile, sep="\t", index=False)
    cannot_gtf.to_csv(ic_ofile, sep="\t", index=False, header=False)
    
    
    



if __name__ == "__main__":
    # parse arguments
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass
    parser = argparse.ArgumentParser(
        description='Create annotation file with the CDS coordinates in several formats',
        epilog='''
Examples:
    python  src/add_cds_coords_appris.py  ...
        ''',
        formatter_class=CustomFormatter )
    parser.add_argument('-ia',  '--input-appris', required=True, help='APPRIS annotation file(s) in GTF format')
    parser.add_argument('-ip',  '--input-cds', required=True, help='Protein table with the CDS coordinates')
    # parser.add_argument('-o',   '--output', required=True, help='Output file with the CDS coordinates in TSV format')
    parser.add_argument('-oc',  '--output-cds', required=True, help='Output file with the CDS coordinates in GTF format')
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

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
    ofile    = args.output
    ic_ofile = args.output_cds
    pg_ofile = args.output_pep
    # ia_ifile = r"S:\U_Proteomica\UNIDAD\Softwares\jmrodriguezc\SANPRO\tmp\202501\zebrafish\appris_method.spade.gtf"
    # ic_ofile = r"S:\U_Proteomica\UNIDAD\Softwares\jmrodriguezc\SANPRO\data\202501\zebrafish\zebrafish_202501.spade.cds.gtf"
    # pg_ofile = r"S:\U_Proteomica\UNIDAD\Softwares\jmrodriguezc\SANPRO\data\202501\zebrafish\zebrafish_202501.spade.pep.gtf"
    # ofile    = r"S:\U_Proteomica\UNIDAD\Softwares\jmrodriguezc\SANPRO\data\202501\zebrafish\zebrafish_202501.spade.tsv"




    logging.info("reading SPADE annotation table...")
    appris_gtf = pd.read_csv(ia_ifile, sep="\t", header=None, low_memory=False, names=GTF_COLUMNS)
    appris_gtf = appris_gtf.astype(str)




    logging.info("extracting SPADE annotation...")
    # extract transcript_id from the attributes column
    attribute_regex = r'.*transcript_id\s"([^"]+)"'
    appris_gtf[['transcript_id']] = appris_gtf['attributes'].str.extract(attribute_regex)
    # extract gene_id from the attributes column
    attribute_regex = r'.*gene_id\s"([^"]+)"'
    appris_gtf[['gene_id']] = appris_gtf['attributes'].str.extract(attribute_regex)
    # extract attributes...
    attribute_regex = r'.*hmm_acc:([^,]+)'
    appris_gtf[['hmm_acc']] = appris_gtf['attributes'].str.extract(attribute_regex)
    # extract attributes...
    attribute_regex = r'.*hmm_name:([^,]+)'
    appris_gtf[['hmm_name']] = appris_gtf['attributes'].str.extract(attribute_regex)
    # extract attributes...
    attribute_regex = r'.*evalue:([^,]+)'
    appris_gtf[['evalue']] = appris_gtf['attributes'].str.extract(attribute_regex)    
    # extract attributes...
    attribute_regex = r'.*pep_start:(\d+)'
    appris_gtf[['pep_start']] = appris_gtf['attributes'].str.extract(attribute_regex)
    # extract attributes...
    attribute_regex = r'.*pep_end:(\d+)'
    appris_gtf[['pep_end']] = appris_gtf['attributes'].str.extract(attribute_regex)
    
    


    logging.info("adding the chr value into CDS coordinate of protein table...")
    # get 'seqname' (chr) and transcript_id without duplicates
    annot_df = appris_gtf.drop_duplicates().rename(columns={
        'seqname':'chr',
        'start':'cds_start',
        'end':'cds_end',
        'strand': 'cds_strand',
        'frame': 'cds_frame',
        'evalue':'pfamscan_evalue'
        })
    # rename feature attribute with the appris annotation
    annot_df['source'] = "SPADE"




    logging.info("create report with CDS coordinates in GTF format...")
    # create the report with the CDS coordinates
    cannot_gtf = annot_df[[
        'chr',
        'source',
        'feature',
        'cds_start',
        'cds_end',
        'score',
        'cds_strand',
        'cds_frame'
    ]].copy()
    # rename
    cannot_gtf = cannot_gtf.rename(columns={
        'cds_strand':'strand',
        'cds_frame':'frame'
        })
    # get attributes
    cannot_gtf['attributes'] = ""
    cc = ['transcript_id','gene_id','pep_start','pep_end','hmm_acc','hmm_name','pfamscan_evalue']
    for c in cc:
        cannot_gtf['attributes'] += c + '=' + annot_df[c] + ';'
    # overwrite values
    cannot_gtf['feature'] = "CDS"
    cannot_gtf['source'] = "GENCODE"




    logging.info("create report with PEPtide coordinates in GTF format...")
    # create the report with the PEP coordinates
    pannot_gtf = annot_df[[
        'transcript_id',
        'source',
        'feature',
        'pep_start',
        'pep_end',
        'score',
        'cds_strand',
        'cds_frame'
    ]].copy()
    # overwrite values
    pannot_gtf['cds_strand'] = "."
    pannot_gtf['cds_frame'] = "."
    # rename
    pannot_gtf = pannot_gtf.rename(columns={
        'transcript_id':'seqname',
        'pep_start':'start',
        'pep_end':'end',
        'cds_strand':'strand',
        'cds_frame':'frame'
        })
    # get attributes
    pannot_gtf['attributes'] = ""
    cc = ['transcript_id','gene_id','chr','cds_start','cds_end','cds_strand','cds_frame','pep_start','pep_end','feature','score','hmm_acc','hmm_name','pfamscan_evalue']
    for c in cc:
        if c == 'score':
            pannot_gtf['attributes'] += 'spade_score=' + annot_df[c] + ';'
        elif c == 'feature':
            pannot_gtf['attributes'] += 'spade_label=' + annot_df[c] + ';'
        else:
            pannot_gtf['attributes'] += c + '=' + annot_df[c] + ';'




    logging.info("create report TSV format...")
    # get columns
    pannot_df = annot_df[[
        'transcript_id',
        'source',
        'feature',
        'pep_start',
        'pep_end',
        'score',
        'cds_strand',
        'cds_frame'
    ]].copy()
    # overwrite values
    pannot_df['cds_strand'] = "."
    pannot_df['cds_frame'] = "."
    # rename
    pannot_df = pannot_df.rename(columns={
        'transcript_id':'seqname',
        'pep_start':'start',
        'pep_end':'end',
        'cds_strand':'strand',
        'cds_frame':'frame'
        })
    # get attributes
    cc = ['transcript_id','gene_id','chr','cds_start','cds_end','cds_strand','cds_frame','pep_start','pep_end','feature','score','hmm_acc','hmm_name','pfamscan_evalue']
    for c in cc:
        if c == 'score':
            pannot_df['spade_score'] = annot_df[c]
        elif c == 'feature':
            pannot_df['spade_label'] = annot_df[c]
        else:
            pannot_df[c] = annot_df[c]





    logging.info("printing output...")
    cannot_gtf.to_csv(ic_ofile, sep="\t", index=False, header=False)
    pannot_gtf.to_csv(pg_ofile, sep="\t", index=False, header=False)
    pannot_df.to_csv(ofile, sep="\t", index=False)
    
    
    



if __name__ == "__main__":
    # parse arguments
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass
    parser = argparse.ArgumentParser(
        description='Create annotation file with the CDS coordinates in several formats',
        epilog='''
Examples:
    python  src/add_cds_coords_spade.py  ...
        ''',
        formatter_class=CustomFormatter )
    parser.add_argument('-ia',  '--input-appris', required=True, help='APPRIS annotation file(s) in GTF format')
    parser.add_argument('-oc',  '--output-cds', required=True, help='Output file with the CDS coordinates in GTF format')
    parser.add_argument('-op',  '--output-pep', required=True, help='Output file with the PEP coordinates in GTF format')
    parser.add_argument('-o',   '--output', required=True, help='Output file with the CDS coordinates in TSV format')
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

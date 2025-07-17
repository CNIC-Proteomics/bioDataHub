# -*- coding: utf-8 -*-
"""
@author: jmrodriguezc
"""

# import global modules
import sys
import argparse
import logging
import pandas as pd



###################
# Parse arguments #
###################

parser = argparse.ArgumentParser(
    description='Add the intersected APPRIS annotations (overlap, PI and NPI CDS) to global annotation file',
    epilog='''Examples:
        
    python  add_intersect_appris.py
      -i  appris.cds_coords.gtf
      -ii appris.cds_coords.intersect.gtf
      -iip appris.cds_coords.unique_PI.gtf
      -iin appris.cds_coords.unique_NPI.gtf
      -oc  appris.cds_coords.overlap.gtf
      -op  appris.pep_coords.overlap.gtf
      -ot  appris.pep_coords.overlap.tsv
    ''',
    formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('-i',    required=True, help='Table with coordinates in GTF format')
parser.add_argument('-ii',   required=True, help='GTF with the intersect coordiantes')
parser.add_argument('-iip',  required=True, help='GTF with the PI coordinates')
parser.add_argument('-iin',  required=True, help='GTF with the NPI coordinates')
parser.add_argument('-oc',   required=True, help='Output file with the CDS coordinates in GTF format')
parser.add_argument('-op',   required=True, help='Output file with the PEP coordinates in GTF format')
parser.add_argument('-ot',   required=True, help='Output file with the PEP coordinates in TSV format')
args = parser.parse_args()

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')


#############
# Constants #
#############
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

APPRIS_COLUMNS = [
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
    'note',
    'gene_name',
    'uniprot_id'
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
    ifile     = args.i
    ifile_ii  = args.ii
    ifile_iip = args.iip
    ifile_iin = args.iin
    ofile_cds = args.oc
    ofile_pep = args.op
    ofile_dat = args.ot
    # ifile     = r"S:\U_Proteomica\UNIDAD\Databases\APPRIS\202501_2\mouse\mouse_202501.appris.cds.gtf"
    # ifile_ii  = r"S:\U_Proteomica\UNIDAD\Databases\APPRIS\202501_2\mouse\mouse_202501.appris.cds_int.gtf"
    # ifile_iip = r"S:\U_Proteomica\UNIDAD\Databases\APPRIS\202501_2\mouse\mouse_202501.appris.cds_u_pi.gtf"
    # ifile_iin = r"S:\U_Proteomica\UNIDAD\Databases\APPRIS\202501_2\mouse\mouse_202501.appris.cds_u_npi.gtf"
    # ofile_cds = r"S:\U_Proteomica\UNIDAD\Databases\APPRIS\202501_2\mouse\mouse_202501.appris.cds_overlap.gtf"
    # ofile_pep = r"S:\U_Proteomica\UNIDAD\Databases\APPRIS\202501_2\mouse\mouse_202501.appris.pep.gtf"
    # ofile_dat = r"S:\U_Proteomica\UNIDAD\Databases\APPRIS\202501_2\mouse\mouse_202501.appris.tsv"


    

    logging.info("reading coordinate table...")
    rep_coords = pd.read_csv(ifile, sep="\t", dtype=str, header=None, low_memory=False)


    

    
    logging.info("reading intersect tables (intesect, PI and NPI)...")
    rep_inter     = pd.read_csv(ifile_ii, sep="\t",  dtype=str, header=None, low_memory=False)
    rep_inter_pi  = pd.read_csv(ifile_iip, sep="\t", dtype=str, header=None, low_memory=False)
    rep_inter_npi = pd.read_csv(ifile_iin, sep="\t", dtype=str, header=None, low_memory=False)

    rep_inter = rep_inter.iloc[:,0:8].drop_duplicates()
    rep_inter_pi = rep_inter_pi.iloc[:,0:8].drop_duplicates()
    rep_inter_npi = rep_inter_npi.iloc[:,0:8].drop_duplicates()




    # # fixing a bug in the case of zebrafish
    # # ValueError: You are trying to merge on object and int64 columns for key '0'. If you wish to proceed you should use pd.concat
    # rep_inter[0] = rep_inter[0].astype('object')
    # rep_inter_pi[0] = rep_inter_pi[0].astype('object')
    # rep_inter_npi[0] = rep_inter_npi[0].astype('object')



    # create empty df that saves the mergining
    rep_merge = pd.DataFrame(columns=['intersect', 'unique_PI', 'unique_NPI'])

    logging.info("add 'intersect' coordinates...")
    # merge based on the given columns add a label when it is in both
    rep_merge['intersect'] = rep_coords.merge(rep_inter, on=[0,1,2,3,4,6,7], how='left', indicator=True)['_merge'] == 'both'
    rep_merge['intersect'] = rep_merge['intersect'].replace({True: 'overlap', False: ''})
    

    logging.info("add 'unique PI' coordinates...")
    # merge based on the given columns add a label when it is in both
    rep_merge['unique_PI'] = rep_coords.merge(rep_inter_pi, on=[0,1,2,3,4,6,7], how='left', indicator=True)['_merge'] == 'both'
    rep_merge['unique_PI'] = rep_merge['unique_PI'].replace({True: 'PI', False: ''})
    
    
    logging.info("add 'unique Non-PI' coordinates...")
    # merge based on the given columns add a label when it is in both
    rep_merge['unique_NPI'] = rep_coords.merge(rep_inter_npi, on=[0,1,2,3,4,6,7], how='left', indicator=True)['_merge'] == 'both'
    rep_merge['unique_NPI'] = rep_merge['unique_NPI'].replace({True: 'NPI', False: ''})

    # concat both dataframes: the report with coordinates and the report with the merged information
    rep_coords = pd.concat([rep_coords, rep_merge], axis=1, join="inner")

    logging.info("add annotations into 'overlap_label' attribute...")
    rep_coords[8] += ";overlap_labels="+rep_coords['intersect']+","+rep_coords['unique_PI']+","+rep_coords['unique_NPI']
    rep_coords[8] = rep_coords[8].str.replace(r';+', ';', regex=True)
    rep_coords[8] = rep_coords[8].str.replace(r',+', ',', regex=True)
    rep_coords[8] = rep_coords[8].str.replace(r',$', '', regex=True)
    rep_coords[8] = rep_coords[8].str.replace(r'=,+', '=', regex=True)

    



    logging.info("printing the CDS output file...")
    rep_coords = rep_coords.iloc[:, 0:9]
    rep_coords.to_csv(ofile_cds, sep="\t", index=False, header=False)
    
    
    
    

    logging.info("obtaining the PEP report from the CDS overlapping...")
    # extract transcript_id from the attributes column
    attribute_regex = r'.*transcript_id=([^;]+);'
    rep_coords['transcript_id'] = rep_coords[8].str.extract(attribute_regex)
    # extract pep_start
    attribute_regex = r'.*pep_start=(\d+)'
    rep_coords['pep_start'] = rep_coords[8].str.extract(attribute_regex)
    # extract pep_end
    attribute_regex = r'.*pep_end=(\d+)'
    rep_coords['pep_end'] = rep_coords[8].str.extract(attribute_regex)




    
    
    logging.info("printing the PEP output file...")
    # create GTF for peptide coordinates
    rep_coords[0] = rep_coords['transcript_id']
    rep_coords[3] = rep_coords['pep_start']
    rep_coords[4] = rep_coords['pep_end']
    rep_coords[6] = '.'
    rep_coords[7] = '.'
    rep_pep_coords = rep_coords.iloc[:, 0:9]
    rep_pep_coords.to_csv(ofile_pep, sep="\t", index=False, header=False)





    logging.info("printing the TSV output file from PEP GTF...")
    # create TSV for peptide coordinates
    rep_pep_coords.columns = GTF_COLUMNS
    # split the attributes column by ";" into multiple columns
    pannot_attrs = rep_pep_coords['attributes'].str.split(';', expand=True)
    # remove labels (like ID=, Parent=, Gene=, Note=) and create columns with their names
    for col in pannot_attrs.columns:
        # Extract the key (e.g., ID, Parent, Gene) and value
        pannot_attrs[col] = pannot_attrs[col].str.strip()
        if "=" in pannot_attrs[col].iloc[0]:
            key = pannot_attrs[col].iloc[0].split('=')[0].strip()
            pannot_attrs[col] = pannot_attrs[col].str.split('=', expand=True)[1]
            pannot_attrs = pannot_attrs.rename(columns={col: key})
    # remove attr column
    rep_pep_coords = rep_pep_coords.drop(columns=['attributes']).reset_index(drop=True)
    pannot_attrs = pannot_attrs.reset_index(drop=True)
    # combine the split attributes back with the original DataFrame
    rep_dat_coords = pd.concat([rep_pep_coords, pannot_attrs], axis=1)
    rep_dat_coords.to_csv(ofile_dat, sep="\t", index=False)
    




if __name__ == "__main__":
    # start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    main(args)
    logging.info('end script')


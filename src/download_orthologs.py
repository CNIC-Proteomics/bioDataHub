#!/usr/bin/python
__author__ = 'jmrodriguezc'
__credits__ = ["Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "2.4"
__maintainer__ = "Jose Rodriguez"
__email__ = "jmrodriguezc@cnic.es"
__status__ = "Development"


import sys
import os
import argparse
import logging
import json
import re
import urllib.request
import urllib.parse

####################
# Global variables #
####################

URL='http://www.ensembl.org/biomart/martservice?query='
QUERY='''
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query  virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="0" count="" datasetConfigVersion="0.6" >
			
	<Dataset name="__SPECIES___gene_ensembl" interface="default" >
		<Filter name="with_hsapiens_homolog" excluded="0"/>
		<Attribute name="ensembl_gene_id" />
		<Attribute name="ensembl_transcript_id" />
		<Attribute name="ensembl_peptide_id" />
		<Attribute name="hsapiens_homolog_ensembl_gene" />
		<Attribute name="hsapiens_homolog_associated_gene_name" />
		<Attribute name="hsapiens_homolog_ensembl_peptide" />
		<Attribute name="hsapiens_homolog_canonical_transcript_protein" />
		<Attribute name="hsapiens_homolog_subtype" />
		<Attribute name="hsapiens_homolog_orthology_type" />
		<Attribute name="hsapiens_homolog_perc_id" />
		<Attribute name="hsapiens_homolog_perc_id_r1" />
		<Attribute name="hsapiens_homolog_goc_score" />
		<Attribute name="hsapiens_homolog_wga_coverage" />
		<Attribute name="hsapiens_homolog_orthology_confidence" />
	</Dataset>
</Query>
'''

#################
# Main function #
#################
def main(args):
    
    logging.info("preparing the query for the given species...")
    # Local folder
    localdir = os.path.dirname(os.path.abspath(__file__))
    # Config file for the species
    with open(os.path.join(localdir, 'config.json')) as f:
       species_cfg = json.load(f)
       
    # check species in the config file
    if args.species in species_cfg.keys() and 'alias' in species_cfg[args.species]:
        # getting the alias
        species = species_cfg[args.species]['alias']

        # preparing the query for the given species...
        query = re.sub('\s*\n\s*', '', QUERY.replace('__SPECIES__',species))
        
        # download orthologs
        logging.info("downloading the human orthologs...")
        q = URL + urllib.parse.quote(query)
        urllib.request.urlretrieve(q, args.outfile)
    
    else:
        logging.error(f"The ({args.species}) species does not exists or does not have an alias")


if __name__ == "__main__":
    # parse arguments
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass
    parser = argparse.ArgumentParser(
        description='Download the human orthologs from Ensembl BioMart',
        epilog='''
Examples:
    python src/download_orthologs.py --species 'pig' --outfile test/biomart_pig.txt
        ''',
        formatter_class=CustomFormatter )
    parser.add_argument('-s',  '--species', required=True, help='First filter based on the species name')
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

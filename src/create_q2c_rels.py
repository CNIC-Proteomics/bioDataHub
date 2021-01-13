#!/usr/bin/python
import sys, argparse, logging
from Bio import SeqIO


__author__ = 'jmrodriguezc'
__credits__ = ["Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.0.2"
__maintainer__ = "Jose Rodriguez"
__email__ = "jmrodriguezc@cnic.es"
__status__ = "Development"


def main(args):
    ''' Main function'''
            
    logging.info("parse the input file with the protein sequences")
    indb = SeqIO.index(args.infasta, "fasta", key_function=lambda rec : rec.split("|")[1])

    print("wait")

if __name__ == "__main__":
    # parse arguments
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass
    parser = argparse.ArgumentParser(
        description='Create the System Biology database from UniProtKB',
        epilog='''
Examples:
    create_q2c_rels.py -f rat_202011_sw.fasta -d rat_202011_sw.categories.tsv
        ''',
        formatter_class=CustomFormatter )
    parser.add_argument('-f',  '--infasta', required=True, help='Protein sequence file in FASTA format')
    parser.add_argument('-d',  '--indb', required=True, help='Tabular file with protein information from UniProt')
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

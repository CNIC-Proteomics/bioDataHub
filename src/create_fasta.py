#!/usr/bin/python
import sys, argparse, logging
import db

__author__ = 'jmrodriguezc'
__credits__ = ["Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "0.0.2"
__maintainer__ = "Jose Rodriguez"
__email__ = "jmrodriguezc@cnic.es"
__status__ = "Development"


def main(args):
    ''' Main function'''
            
    logging.info("create db_creator object")
    w = db.creator(args.species, args.outfile, args.filt)

    logging.info("download fasta file")
    w.download_fasta_dbs(args.filt, args.rem_dup)



if __name__ == "__main__":
    # parse arguments
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass
    parser = argparse.ArgumentParser(
        description='Dowloads FASTA sequence file determining the source of database, UniProtKB or the UniProt Proteome',
        epilog='''
Examples:
    create_fasta.py -s pig   -f uni-sw     -o test
    create_fasta.py -s mouse -f pro-sw-tr  -o test
        ''',
        formatter_class=CustomFormatter )
    parser.add_argument('-s',  '--species', required=True, help='First filter based on the species name')
    parser.add_argument('-o',  '--outfile', required=True, help='Output file')
    parser.add_argument('-f',  '--filt',  required=True, choices=["pro-sw-tr","pro-sw","uni-sw-tr","uni-sw"], help='Determines which kind of database you download')
    parser.add_argument('-d',  '--rem_dup', default=False, action='store_true', help="Remove duplicated sequences")
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

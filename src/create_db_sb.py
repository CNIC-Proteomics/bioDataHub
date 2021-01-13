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
    w = db.creator(args.species, args.outdir, args.filt, args.rem_dup)

    logging.info("download raw file")
    w.download_raw_dbs(args.filt)

    logging.info("create qreport")
    output = w.create_qreport()

    logging.info('print database file')
    w.to_file(output)


if __name__ == "__main__":
    # parse arguments
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass
    parser = argparse.ArgumentParser(
        description='Create the System Biology database from UniProtKB',
        epilog='''
Examples:
    create_db_sb.py -s pig -o test
    create_db_sb.py -s mouse -f -t gene -o test
        ''',
        formatter_class=CustomFormatter )
    parser.add_argument('-s',  '--species', required=True, help='First filter based on the species name')
    parser.add_argument('-o',  '--outdir', required=True, help='Directory where the database will be saved')
    parser.add_argument('-f',  '--filt',  default="sw-tr", choices=["sw-tr","sw","tr"], help='Directory where the database will be saved')    
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

#!/usr/bin/python

__author__ = 'jmrodriguezc'
__credits__ = ["Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "2.4"
__maintainer__ = "Jose Rodriguez"
__email__ = "jmrodriguezc@cnic.es"
__status__ = "Development"


import sys
import argparse
import logging
from Bio import SeqIO
from difflib import SequenceMatcher


def similar(a, b):
    return SequenceMatcher(None, a, b).ratio()

def sequence_cleaner(sequences, fasta_file, por_n=100):
    # Using the Biopython fasta parse we can read our fasta input
    for seq_record in SeqIO.parse(fasta_file, "fasta"):
        # Take the current sequence
        sequence = str(seq_record.seq).upper()
        # Check if the current sequence is according to the user parameters
        # If the sequence passed in the test "is it clean?" and it isn't in the
        # hash table, the sequence and its id are going to be in the hash
        if sequence not in sequences:
            sequences[sequence] = [seq_record.id]
        # If it is already in the hash table, we're just gonna concatenate the ID
        # of the current sequence to another one that is already in the hash table
        else:
            sequences[sequence].append(seq_record.id)

def main(args):

    # Create our hash table to add the sequences
    sequences = {}    
    [ sequence_cleaner(sequences,ifile) for ifile in args.infiles ]

    # Create a file in the same directory where you ran this script
    with open(args.outfile, "w+") as output_file:
        # Just read the hash table and write on the file as a fasta format
        for sequence in sequences:
            output_file.write(">" + "__".join(sequences[sequence]) + "\n" + sequence + "\n")
            

if __name__ == "__main__":
    # parse arguments
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass
    parser = argparse.ArgumentParser(
        description='Remove redundance of FASTA files',
        epilog='''
Examples:
    remove_redundances.py -ii 
        ''',
        formatter_class=CustomFormatter )
    parser.add_argument('-i',  '--infiles', required=True, nargs='+', help='FASTA files')
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

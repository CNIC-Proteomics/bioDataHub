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
import gzip
import shutil
import requests
import urllib.request
import urllib.parse

####################
# Global variables #
####################

URL='https://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/'
# URL='https://apprisws.bioinfo.cnio.es/pub/releases/2024_06.v48/datafiles/'


###################
# Local functions #
###################
def download_appris_dbs(url, exts, odir):
    # Send an HTTP GET request to the server directory
    response = requests.get(url)
    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        # Retrieves the hrefs from the HTML content of the page
        hrefs = re.findall(r'href=\"([^\"]*)\"', response.text)
        # Check if the filename ends with one of the target extensions
        hrefs = [ h for h in hrefs for e in exts if h.endswith(e) ]
        # Loop through the hrefs
        for href in hrefs:
            try:
                # Construct the full URL for the file
                href_url = url + '/'+ href
                # Download the file
                db_dat = odir +'/'+ os.path.basename(href)
                logging.info("get "+href_url+" > "+db_dat)
                urllib.request.urlretrieve(href_url, db_dat)
            except Exception as exc:
                logging.warning(f"failed to dowload {href}: {exc}")
            try:
                if href.endswith('.gz'):
                    db_dat2 = os.path.join(odir, '.'.join(os.path.basename(href).split('.')[:-1]))
                    logging.info("uncompressing "+db_dat2)
                    with gzip.open(db_dat, 'rb') as f_in:
                        with open(db_dat2, 'wb') as f_out:
                            shutil.copyfileobj(f_in, f_out)
            except Exception as exc:
                logging.warning(f"failed uncompressing the file {db_dat}: {exc}")                    

#################
# Main function #
#################
def main(args):
    
    logging.info("preparing the workspace...")
    os.makedirs(args.outdir, exist_ok=False)

    logging.info("preparing the query for the given species...")
    # Local folder
    localdir = os.path.dirname(os.path.abspath(__file__))
    # Config file for the species
    with open(os.path.join(localdir, 'config.json')) as f:
       species_cfg = json.load(f)
       
    # check species in the config file
    if args.species in species_cfg.keys() and 'assembly' in species_cfg[args.species]:
        logging.info("preparing the base url...")
        # getting species info
        assembly = species_cfg[args.species]['assembly']
        scientific = species_cfg[args.species]['scientific']
        # preparing base url
        name = scientific.lower().replace(' ','_')
        baseurl = URL +name+'/'+assembly
        
        logging.info("downloading the method annotations...")
        download_appris_dbs(baseurl, [
            '.crash.gtf.gz',
            '.firestar.gtf.gz',
            '.matador3d.gtf.gz',
            '.spade.gtf.gz',
            '.thump.gtf.gz',
            '.appris.gtf.gz',
            '.pannot.gtf'],
            args.outdir)
    
    else:
        logging.error(f"The ({args.species}) species does not exists or does not have an assembly")


if __name__ == "__main__":
    # parse arguments
    class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter, argparse.RawTextHelpFormatter):
        pass
    parser = argparse.ArgumentParser(
        description='Downloads the annotations for the APPRIS methods that locate the annotation in a specific region of the protein',
        epilog='''
Examples:
    python src/download_appris.py --species 'pig' --outdir test/appris
        ''',
        formatter_class=CustomFormatter )
    parser.add_argument('-s',  '--species', required=True, help='First filter based on the species name')
    parser.add_argument('-o',  '--outdir', required=True, help='Output folder')
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

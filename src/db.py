import sys, os, logging
import urllib.request
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry
import datetime
import re
from concurrent.futures import ThreadPoolExecutor
import gzip
import shutil
import json
import pandas as pd
import numpy as np
import time
import socket
from Bio import SwissProt
from Bio import SeqIO
from Bio.KEGG import REST
from OpenSSL import SSL, crypto


class creator:
    # Local folder
    LOCAL_DIR = os.path.dirname(os.path.abspath(__file__))

    # Config file for the species
    with open(os.path.join(LOCAL_DIR, 'config.json')) as f:
       SPECIES_CFG = json.load(f)

    # https://rest.uniprot.org/uniprotkb/stream?format=fasta&includeIsoform=true&query=%28taxonomy_id%3A10090%29%20AND%20%28reviewed%3Atrue%29
    # https://rest.uniprot.org/uniprotkb/stream?format=fasta&includeIsoform=true&query=%28%28proteome%3AUP000000589%29%29
    URL_UNIPROT = 'https://rest.uniprot.org/uniprotkb/stream?'
    URL_UNIPROT += 'includeIsoform=true&' # include all isoforms    
    # URL_CORUM   = 'http://mips.helmholtz-muenchen.de/corum/download/allComplexes.json.zip' #It doesn't work :-(
    # URL_CORUM   = os.path.join(LOCAL_DIR, '../cached/allComplexes.json.zip')
    URL_CORUM   = 'https://mips.helmholtz-muenchen.de/fastapi-corum/public/file/download_current_file?file_id=complete&file_format=json'
    URL_PANTHER = 'http://data.pantherdb.org/ftp/sequence_classifications/current_release/PANTHER_Sequence_Classification_files/'
    URL_APPRIS  = 'https://apprisws.bioinfo.cnio.es/pub/current_release/datafiles/'
    cRAP_FILE   = os.path.join(LOCAL_DIR, '../cached/cRAP/crap.modified.fasta')
    # Column name with the cross-reference id
    XID = 'Protein'
    # Meta terms of Protein
    META = ['xref_UniProt_Name','Protein','Gene','Species','Length','Description','Comment_Line','prot_UniProt_Class']
    # Xreferences terms of Protein
    XTERMS = [
        ('Ensembl',  ['xref_Ensembl_protId','xref_Ensembl_transcId','xref_Ensembl_GeneId','Protein']),
        ('RefSeq',   ['xref_RefSeq_protId','xref_RefSeq_transcId','Protein']),
        ('CCDS',     ['xref_CCDS','Protein'])
    ]
    # Category terms of Protein
    CTERMS = [
        ('GO',       [('cat_GO_C','C'),('cat_GO_F','F'),('cat_GO_P','P')]),
        ('KEGG',     [('cat_KEGG','')]),
        ('PANTHER',  [('cat_PANTHER','')]),
        ('Reactome', [('cat_Reactome','')]),
        ('CORUM',    [('cat_CORUM','')]),
        ('MIM',      [('cat_OMIM','')]),
        ('DrugBank', [('cat_DrugBank','')])
    ]
    # APPRIS terms of Protein
    ATERMS = [
        ('APPRIS',   [('APPRIS Annotation','')]),
        ('TRIFID',   [('norm_trifid_score','')]),
        ('CORSAIR',   [('corsair_score','')])
    ]
    HEADERS = [ h for h in META] + [ h for i in XTERMS for h in i[1][:-1] ] + [ h[0] for i in CTERMS for h in i[1] ] + [ h[0] for i in ATERMS for h in i[1] ]
    TIME = datetime.datetime.now().strftime("%Y%m")


    '''
    Creates the databases
    '''
    def __init__(self, s, o, f=None):
        # assign species
        if s in self.SPECIES_CFG.keys():
            self.species = s
            self.proteome_id = self.SPECIES_CFG[self.species]['proteome']
            self.kegg_id = self.SPECIES_CFG[self.species]['kegg_organism']
            self.scientific = self.SPECIES_CFG[self.species]['scientific']
            self.taxonomy = self.SPECIES_CFG[self.species]['taxonomy']
            self.assembly = self.SPECIES_CFG[self.species]['assembly'] if 'assembly'in self.SPECIES_CFG[self.species] else None
        else:
            sys.exit( "ERROR: Species parameter has been not found. Try with: "+", ".join(self.SPECIES_CFG.keys()) )
        
        # prepare temporal file
        self.TMP_DIR = self.LOCAL_DIR +'/../tmp/'+ self.TIME +'/'+ self.species
        os.makedirs(self.TMP_DIR, exist_ok=True)
        logging.debug(f"TMP_DIR: {self.TMP_DIR}")

        # prepare the output directory if does not exist
        self.outdir = os.path.dirname(o)
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir, exist_ok=True)
        
        # define output file
        self.outfile = o

        # define output files for "create_sb"
        self.db_uniprot = self.TMP_DIR +'/uniprot.dat'
        # self.db_uniprot = self.TMP_DIR +'/test.dat'
        # self.db_uniprot = self.TMP_DIR +'/../../../test/test_1033.dat'
        
        self.db_fasta = self.TMP_DIR +'/proteins.fasta'
        self.db_corum   = self.TMP_DIR +'/'+ ".".join(['allComplexes', 'json'])
        self.db_panther = self.TMP_DIR +'/panther.dat'
        self.db_kegg    = self.TMP_DIR +'/kegg.dat'
        self.db_appris  = self.TMP_DIR +'/appris.dat'
        self.db_trifid  = self.TMP_DIR +'/trifid.dat'
        self.db_corsair  = self.TMP_DIR +'/corsair.dat'
        self.CACHE_DIR = os.path.join(self.TMP_DIR, "kegg_cache")
        os.makedirs(self.CACHE_DIR, exist_ok=True)


    def _download_file(self, url, dest_path, retries=5, backoff_factor=0.3, chunk_size=16*1024):
        '''
        Download file from URL with retries
        '''
        try:
            session = requests.Session()
            retry = Retry(
                total=retries,
                read=retries,
                connect=retries,
                backoff_factor=backoff_factor,
                status_forcelist=(500, 502, 504),
            )
            adapter = HTTPAdapter(max_retries=retry)
            session.mount('http://', adapter)
            session.mount('https://', adapter)

            with session.get(url, headers={"Accept-Encoding": "identity"}, stream=False) as r: # disables gzip/deflate compression
                r.raise_for_status()
                with open(dest_path, 'wb') as f:
                    for chunk in r.iter_content(chunk_size=chunk_size, decode_unicode=False):
                        if chunk:  # filter out keep-alive new chunks
                            f.write(chunk)
        except requests.exceptions.RequestException as e:
            logging.error(f"Request error: {e}")
            raise
                
    def _delete_tmp_dir(self, dir):
        files = [ f for f in os.listdir(dir) ]
        for f in files:
            try:
                os.remove(os.path.join(dir, f))
            except Exception as e:
                logging.error(e)

    def _add_crap(self, infile, outfile=None):
        '''
        Add the cRAP database
        '''
        infile2 = self.cRAP_FILE
        tmp_outfile = outfile+'.txt'
        with open(tmp_outfile, 'w') as fo, open(infile, 'r') as f1, open(infile2, 'r') as f2:
            fo.write("".join("{}{}".format(f1.read(),f2.read()) ))
        # remove obsolete output file
        if os.path.isfile(outfile):
            os.remove(outfile)
        # rename the temporal file
        os.rename(tmp_outfile, outfile)
        return None

    def _remove_duplicates(self, infile, outfile=None):
        '''
        Remove duplicated sequences
        '''
        seqs = dict()
        precords = SeqIO.parse(infile, "fasta")
        records = SeqIO.to_dict(precords)
        # read all sequences
        # delete the trEMBL sequences based on the sorted id's (tr prefix)
        for i in list(records):
            record = records[i]
            s = str(record.seq)
            if s in seqs:
                l = [seqs[s], i]
                # check if there is a SwissProt protein and TrEMBL protein
                p = ['sp|','tr|']
                if all(x in " ".join(l) for x in p):
                    logging.warning( "duplicated sequences: {}".format(",".join(l)) , exc_info=False)
                    # remove all trEMBL proteins that are duplicated. remove multiple keys from dictionary
                    logging.warning( "deleting the trEMBL sequences", exc_info=False)
                    [ records.pop(i) for i in l if 'tr|' in i ]
            else:
                seqs[s] = i
        # write file if apply
        if outfile:
            with open(outfile, 'w') as handle:
                SeqIO.write(records.values(), handle, 'fasta')
        return None

    def _fetch_kegg_batch(self, entries):
            # check if all entries cached
            uncached = []
            batch_cache_files = {}
            for e in entries:
                f = os.path.join(self.CACHE_DIR, e.replace(":", "_") + ".txt")
                batch_cache_files[e] = f
                if not os.path.exists(f):
                    uncached.append(e)

            # if all cached, load from disk
            if len(uncached) == 0:
                data = {}
                for e, f in batch_cache_files.items():
                    with open(f) as fh:
                        data[e] = fh.read()
                return data

            # otherwise download entries
            query = "+".join(uncached)
            time.sleep(0.25)
            try:
                raw = REST.kegg_get(query).read()
            except Exception as exc:
                print(f"[ERROR] Failed batch {uncached}: {exc}")
                return {e: None for e in entries}
            blocks = raw.strip().split("///")
            result = {}
            for block in blocks:
                block = block.strip()
                if not block:
                    continue
                first_line = block.split("\n", 1)[0]
                entry_id = self.kegg_id + ":" + first_line.split()[1]
                f = batch_cache_files[entry_id]
                with open(f, "w") as fh:
                    fh.write(block + "\n///\n")
                result[entry_id] = block + "\n///\n"
            # mark missing entries
            for e in entries:
                if e not in result:
                    result[e] = None
            return result

    def _download_kegg(self, db):
        result = {}
        batches = [db[i:i+10] for i in range(0, len(db), 10)] # 10 is max recommended batch size
        with ThreadPoolExecutor(max_workers=6) as ex:
            for i, batch_entries in enumerate(batches):
                if i+1 == 1 or i+1 == len(batches) or (i+1) % 10 == 0:
                    logging.debug(f"Downloading KEGG batch {i+1}/{len(batches)}")
                batch_result = ex.submit(self._fetch_kegg_batch, batch_entries).result()
                result.update(batch_result)
        return result

    def download_fasta_dbs(self, filt=None, d=None):
        '''
        Download the fasta file
        '''        
        # filter by SwissProt (Reviwed) if apply
        if not os.path.isfile(self.outfile):
            # create the query
            query = ''
            if filt and filt.startswith("pro"): # filter by proteome
                query = 'proteome:'+ self.proteome_id
            else: # by default filter by organism
                query = 'taxonomy_id:'+self.taxonomy
            if filt and filt.endswith("sw"): # filter by SwissProt
                query += '%20AND%20reviewed=true'
            url = self.URL_UNIPROT + f"query=({query})"
            url += '&format=fasta'
            logging.info("get "+url+" > "+self.outfile)
            # urllib.request.urlretrieve(url, self.outfile)
            self._download_file(url, self.outfile)
        else:
            logging.info('cached uniprot')
            
        # remove duplicate sequences
        if d:
            logging.info("remove duplicate sequences")
            self._remove_duplicates(self.outfile, self.outfile)
            
        # DEPRECATED
        # # adding the cRAP
        # logging.info("adding the cRAP database")
        # self._add_crap(self.outfile, self.outfile)
            
    
    def _download_appris_db(self, url, name, ofile):
        # Send an HTTP GET request to the server directory
        response = requests.get(url)
        # Check if the request was successful (status code 200)
        if response.status_code == 200:
            # Retrieves the hrefs from the HTML content of the page
            hrefs = re.findall(r'href=\"([^\"]*)\"', response.text)
            # Check if the filename ends with one of the target extensions
            href = [ h for h in hrefs if h == name ]
            if len(href) == 1:
                href = href[0]
                # rename the output if it is compressed file
                if href.endswith('.gz'):
                    ofile = f"{ofile}.gz"
                # dowload the file
                try:
                    # Construct the full URL for the file
                    href_url = url + '/'+ href
                    # Download the file
                    logging.info("get "+href_url+" > "+ofile)
                    # urllib.request.urlretrieve(href_url, ofile)
                    self._download_file(href_url, ofile)
                except Exception as exc:
                    logging.warning(f"failed to dowload {href}: {exc}")
                # uncompress if applied
                try:
                    if href.endswith('.gz'):
                        ofile2 = ofile.replace('.gz','')
                        logging.info("uncompressing "+ofile2)
                        with gzip.open(ofile, 'rb') as f_in:
                            with open(ofile2, 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)
                except Exception as exc:
                    logging.warning(f"failed uncompressing the file {ofile}: {exc}")
            else:
                logging.error(f"failed dowloading {name}")
        else:
            logging.error(f"failed dowloading {name}")

    def getPEMFile(self, host):
        '''
        Get the certificate to download CORUM
        '''

        dst = (host, 443)
        ctx = SSL.Context(SSL.TLS_CLIENT_METHOD)
        s = socket.create_connection(dst)
        s = SSL.Connection(ctx, s)
        s.set_connect_state()
        s.set_tlsext_host_name(str.encode(dst[0]))
        s.sendall(str.encode('HEAD / HTTP/1.0\n\n'))
        peerCertChain = s.get_peer_cert_chain()
        pemFile = ''

        for cert in peerCertChain:
            pemFile += crypto.dump_certificate(crypto.FILETYPE_PEM, cert).decode("utf-8")

        return pemFile

    def download_raw_dbs(self, filt=None):
        '''
        Download the raw databases
        '''
        
        # UniProt Fasta
        if not os.path.isfile(self.db_fasta):
            url = self.URL_UNIPROT +'query=(taxonomy_id:'+self.taxonomy+')'
            url += '&format=fasta'
            logging.info("get "+url+" > "+self.db_fasta)
            # urllib.request.urlretrieve(url, self.db_fasta)
            self._download_file(url, self.db_fasta)
        else:
            logging.info('cached uniprot fasta')
        logging.info('create report from fasta file with the UniProt accession')
        with open(self.db_fasta, 'r') as f:
            fasta_txt = f.read()
        # create a dictionary with the protein accession and the value is a dict(comment, protein_length)
        b = [ a[1:].split('\n') if a.startswith('>') else a.split('\n') for a in fasta_txt.split('\n>') ]
        self.db_fasta_seqio = dict([ (c[0].split("|")[1], {'comm': '>'+c[0], 'len': len("".join(c[1:])) }) for c in b if len(c) > 1])

        # UniProt Data
        if not os.path.isfile(self.db_uniprot):
            url = self.URL_UNIPROT +'query=(taxonomy_id:'+self.taxonomy+')'
            url += '&format=txt'
            logging.info("get "+url+" > "+self.db_uniprot)
            # urllib.request.urlretrieve(url, self.db_uniprot)
            self._download_file(url, self.db_uniprot)
        else:
            logging.info('cached uniprot data')
        
        # CORUM
        # download all complexes file (using the same name)
        if not os.path.isfile(self.db_corum):
            url = self.URL_CORUM
            db_dat = self.TMP_DIR +'/'+ 'allComplexes.json'
            logging.debug("get "+url)
            cert_corum = self.getPEMFile(self.URL_CORUM.split("/")[2])
            cert_file = self.TMP_DIR +'/'+ 'corum_certificate.pem'
            with open(cert_file, "w") as f:
                f.write(cert_corum)
            response = requests.get(url, verify=cert_file)
            response.raise_for_status()          
            with open(db_dat, "wb") as f:
                f.write(response.content)
        else:
            logging.info('cached corum')
        
        # PANTHER
        # get the list of species and extract the file name
        if not os.path.isfile(self.db_panther):
            url = self.URL_PANTHER
            try:
                result = urllib.request.urlopen(url).read().decode('utf-8')
                if result:
                    pattern = re.search(r'\s*(PTHR[^\_]*\_'+self.species+')', result, re.I | re.M)
                    if pattern:
                        url = self.URL_PANTHER + pattern[1]
                        logging.info("get "+url+" > "+self.db_panther)
                        # urllib.request.urlretrieve(url, self.db_panther)
                        self._download_file(url, self.db_panther)
            except Exception as exc:
                logging.warning(f"panther url does not exist: {exc}")
        else:
            logging.info('cached panther')

        # KEGG
        # download in batches
        if not os.path.isfile(self.db_kegg):
            db = REST.kegg_list(database=self.kegg_id).read()
            db = [i.split("\t")[0] for i in db.split("\n")]
            db = [i for i in db if i.startswith(self.kegg_id)]
            result = self._download_kegg(db)
            with open(self.db_kegg, "w") as out:
                for e in db:
                    if result[e]:
                        out.write(result[e])
                    else:
                        logging.warning(f"kegg - failed to retrieve: {e}\n")
        else:
            logging.info('cached kegg')
        shutil.rmtree(self.CACHE_DIR)
    
        # APPRIS
        # get the list of species and extract the file name
        if not os.path.isfile(self.db_appris):
            if self.assembly:
                name = self.scientific.lower().replace(' ','_')
                url = self.URL_APPRIS +name+'/'+self.assembly
                try:
                    self._download_appris_db(url, 'appris_data.principal.txt', self.db_appris)
                except Exception as exc:
                    logging.warning(f"failed dowloading the files: {exc}")
            else:
                logging.warning("assembly does not exist for appris annotation")
        else:
            logging.info('cached appris')
            

        # TRIFID
        # get the list of species and extract the file name
        if not os.path.isfile(self.db_trifid):
            if self.assembly:
                name = self.scientific.lower().replace(' ','_')
                url = self.URL_APPRIS +name+'/'+self.assembly
                try:
                    self._download_appris_db(url, 'appris_method.trifid.txt', self.db_trifid)
                except Exception as exc:
                    logging.warning(f"failed dowloading the files: {exc}")
            else:
                logging.warning("assembly does not exist for trifid annotation")
        else:
            logging.info('cached trifid')
            
            
        # CORSAIR
        # get the list of species and extract the file name
        if not os.path.isfile(self.db_corsair):
            if self.assembly:
                name = self.scientific.lower().replace(' ','_')
                url = self.URL_APPRIS +name+'/'+self.assembly
                try:
                    self._download_appris_db(url, 'appris_method.corsair.gtf.gz', self.db_corsair)
                except Exception as exc:
                    logging.warning(f"failed dowloading the files: {exc}")
            else:
                logging.warning("assembly does not exist for corsair annotation")
        else:
            logging.info('cached corsair')


    
    def create_qreport(self):
        '''
        Create protein report
        '''        
        if self.db_uniprot:
            # create reports from external data ---
            logging.info('create reports from external data...')
            corum_json = None
            panther_df = None
            kegg_dict = None
            appris_df = pd.DataFrame()
            trifid_df = pd.DataFrame()
            corsair_df = pd.DataFrame()
            
            if os.path.isfile(self.db_corum):
                with open(self.db_corum, 'r') as f:
                    corum_json = json.load(f)
                    
            if os.path.isfile(self.db_panther):
                panther_df = pd.read_csv(self.db_panther, sep="\t", dtype=str, header=None, low_memory=False)

            if os.path.isfile(self.db_kegg):
                with open(self.db_kegg, 'r') as f:
                    kegg_dict = f.read()
                kegg_dict = [i for i in kegg_dict.split("///")]
                kegg_dict = {self.kegg_id + ":" + i.split()[1]: i for i in kegg_dict[:-1]}
                
            if os.path.isfile(self.db_appris):
                try:                    
                    appris_df = pd.read_csv(self.db_appris, sep="\t", dtype=str, usecols=['Gene ID','Transcript ID','APPRIS Annotation','MANE'], low_memory=False)
                except Exception as exc:
                    logging.warning(f"reading the appris database: {exc}")
                    # could be that the columns do not exist
                    try:                    
                        appris_df = pd.read_csv(self.db_appris, sep="\t", dtype=str, header=None, low_memory=False)
                        appris_df.columns = ['Gene name', 'Gene ID','Transcript ID','CCDS ID','APPRIS Annotation']
                        appris_df = appris_df[['Gene ID','Transcript ID','APPRIS Annotation']]
                    except Exception as exc:
                        # could be that the columns do not exist
                        logging.error(f"reading the appris database: {exc}")
                        
            if os.path.isfile(self.db_trifid):
                trifid_df = pd.read_csv(self.db_trifid, sep="\t", dtype=str, usecols=['gene_id','transcript_id','norm_trifid_score'], low_memory=False)
                
            if os.path.isfile(self.db_corsair):
                df = pd.read_csv(self.db_corsair, sep="\t", dtype=str, header=None, low_memory=False)
                # create a column with the gene_id and transcriot_id
                df_notes = pd.DataFrame()                
                col_attr = df.columns[-1] # get the name of last column (attributes)
                df_notes['ensembl_gene_id'] = df[col_attr].str.extract(r'gene_id \"([^\"]*)\"', expand=False)
                df_notes['ensembl_transc_id']= df[col_attr].str.extract(r'transcript_id \"([^\"]*)\"', expand=False)
                # join with the rest of information
                df = pd.concat([df,df_notes], axis=1)
                df = df.rename(columns={5:'corsair_score'})
                # extract the id columns and the corsair score (from GTF format)
                corsair_df = df[['ensembl_gene_id','ensembl_transc_id','corsair_score']]
                
            
            logging.info("start the extraction of data...")
            ddf = []
            for rec in SwissProt.parse(open(self.db_uniprot)):
                ddf.append( self._create_qreport(rec,corum_json,panther_df,appris_df,trifid_df,corsair_df) )
            ddf = pd.concat(ddf, axis=0)

        return ddf


    def _create_qreport(self, record, corum_json, panther_df, appris_df, trifid_df, corsair_df):
        '''
        Create protein report
        '''
        # Efficient way to unnest (explode) multiple list columns in a pandas DataFrame
        def _explode(df, lst_cols, fill_value=''):
            # make sure `lst_cols` is a list
            if lst_cols and not isinstance(lst_cols, list):
                lst_cols = [lst_cols]
            # all columns except `lst_cols`
            idx_cols = df.columns.difference(lst_cols)
        
            # calculate lengths of lists
            lens = df[lst_cols[0]].str.len()
        
            if (lens > 0).all():
                # ALL lists in cells aren't empty
                return pd.DataFrame({
                    col:np.repeat(df[col].values, df[lst_cols[0]].str.len())
                    for col in idx_cols
                }).assign(**{col:np.concatenate(df[col].values) for col in lst_cols}) \
                  .loc[:, df.columns]
            else:
                # at least one list in cells is empty
                return pd.DataFrame({
                    col:np.repeat(df[col].values, df[lst_cols[0]].str.len())
                    for col in idx_cols
                }).assign(**{col:np.concatenate(df[col].values) for col in lst_cols}) \
                  .append(df.loc[lens==0, idx_cols]).fillna(fill_value) \
                  .loc[:, df.columns]
          
        
        # extract main info ---
        name = record.entry_name
        acc = record.accessions[0]
        # pattern = re.search(r'Name=([^\s|\;]*)', record.gene_name, re.I | re.M)
        # gene = pattern[1] if pattern else record.gene_name
        gene = "".join([ re.sub(r'\s+.*$', '', r['Name']) if 'Name' in r else re.sub(r'\s+.*$', '', r['ORFNames'][0]) if 'ORFNames' in r else '' for r in record.gene_name ])
        pattern = re.search(r'[RecName|SubName]: Full=([^\;|\{]*)', record.description, re.I | re.M)
        dsc = pattern[1] if pattern else record.description
        dclass = record.data_class
        pattern = re.search(r'([\w|\s]*)\s+\(\w*\)', record.organism, re.I | re.M)
        species = pattern[1] if pattern else record.organism
        # extract isoforms IDs and the displayed isoform
        altprod = [c for c in record.comments if 'ALTERNATIVE PRODUCTS:' in c]
        IsoIds = [acc]
        IsoDisplayed = None
        if altprod:
            ap = re.split(r'\s*;\s*', altprod[0])
            x = [ ( a.replace('IsoId=','') , ap[i+1].replace('Sequence=','') ) for i,a in enumerate(ap) if a.startswith('IsoId=')] # list of tuples (isoId=,Sequence=)
            x = [ ( re.match(r'[^,]*',y[0])[0], re.match(r'[^,]*',y[1])[0] ) for y in x ] # if multiple Iso ids, get the fisrt id until comma
            IsoIds += [ i[0] for i in x if i[1].startswith('VSP_') ] # get the list of
            z =  [ i[0] for i in x if i[1].startswith('Displayed') ]
            if z: IsoDisplayed = z[0]
        # for each isoform id
        # get the comment and length of isoform
        comms = []
        lengths = []
        for i in IsoIds:
            # get the comm and len from iso_id
            c = self.db_fasta_seqio[i]['comm'] if i in self.db_fasta_seqio else ''
            comms.append(c)
            l = self.db_fasta_seqio[i]['len'] if i in self.db_fasta_seqio else ''
            lengths.append(l)        
            
        
        # create a dataframe with the Metadata information ---
        # UniProt accesion isoform is the index
        # ['xref_UniProt_Name','Protein','Gene','Species','Length','Description','Comment_Line','prot_UniProt_Class']
        df1 = pd.DataFrame(columns=self.META, data=[[name,IsoIds,gene,species,lengths,dsc,comms,dclass]])
        df1 = _explode(df1, lst_cols=['Protein','Length','Comment_Line'])
        df1.set_index(self.XID, inplace=True)
        
        # create cross-references data ---
        # filter by given list of terms
        # create dictionary with all common Xreferences
        terms_dbs = [i[0] for i in self.XTERMS]
        rs = [x for x in record.cross_references if x[0] in terms_dbs ]
        rcross = dict()
        for r in rs:
            rcross.setdefault(r[0], []).append(r[1:])
        
        
        # extract the xreference data ---
        for terms in self.XTERMS:
            xdb = terms[0]
            xcols = terms[1]
            xvals = []
            if xdb in rcross:
                rconts = rcross[xdb]
                if   xdb == "Ensembl":
                    for rcont in rconts:
                        xp = rcont[1]
                        xt = rcont[0]
                        # extract the isoform id from the last tuple (gene+isoform)
                        y = re.findall(r"\[([^\]]*)\]", rcont[2])
                        i = acc if not y or y[0] == IsoDisplayed else y[0]
                        # remove the isoform id from the last tuple
                        xg = re.sub(r"\.\s+.*$",'',rcont[2])
                        # remote the version in the identifiers
                        xp = re.sub(r"\.\d+$",'', xp)
                        xt = re.sub(r"\.\d+$",'', xt)
                        xg = re.sub(r"\.\d+$",'', xg)
                        # append
                        xvals.append([xp,xt,xg,i])                            
                elif xdb == "RefSeq":
                    for rcont in rconts:
                        xp = rcont[0]
                        # extract the isoform id from the last tuple (transc+isoform)
                        y = re.findall(r"\[([^\]]*)\]", rcont[1])
                        i = acc if not y or y[0] == IsoDisplayed else y[0]
                        # remove the isoform id from the last tuple
                        xt = re.sub(r"\.\s+.*$",'',rcont[1])
                        # remote the version in the identifiers
                        xp = re.sub(r"\.\d+$",'', xp)
                        xt = re.sub(r"\.\d+$",'', xt)
                        # append
                        xvals.append([xp,xt,i])
                elif xdb == "CCDS":
                    for rcont in rconts:
                        xp = rcont[0]
                        # extract the isoform id from the last tuple
                        y = re.findall(r"\[([^\]]*)\]", rcont[1])
                        i = acc if not y or y[0] == IsoDisplayed else y[0]
                        xvals.append([xp,i])
            else:
                # create empty dict with the name of colum
                xvals = [[np.nan] for x in xcols]
                xvals = list(map(list, zip(*xvals)))
            # create dataframe with the Xreferece information
            df2 = pd.DataFrame(columns=xcols, data=xvals)
            if not df2.dropna(how='all').empty:
                # check if UniProt accession does not exit
                # we add the given Isoforms ids
                if df2[self.XID].dropna(how='all').empty:
                    df2[self.XID] = IsoIds                        
                df2.set_index(self.XID, inplace=True)
                # join using the index which is the UniProt accession of isoform
                df1 = df1.join(df2, how='outer')
        
        
        # create cross-references data ---
        # filter by given list of terms
        # create dictionary with all common Xreferences
        terms_dbs = [i[0] for i in self.CTERMS]
        rs = [x for x in record.cross_references if x[0] in terms_dbs ]
        rcross = dict()
        for r in rs:
            rcross.setdefault(r[0], []).append(r[1:])

        # extract the category data from UniProt ---
        for terms in self.CTERMS:
            xdb = terms[0]
            xpats = terms[1]
            xcols,xvals = [],[]
            if xdb in rcross:
                rconts = rcross[xdb]
                if xdb == "GO":
                    (xcols, xvals) = self._extract_cat_go(rconts, xpats)
                elif xdb == "KEGG": # remote access
                    (xcols, xvals) = self._extract_cat_kegg(kegg_dict, rconts, xpats)
                elif xdb == "PANTHER":
                    (xcols, xvals) = self._extract_cat_panther(panther_df, rconts, xpats, acc)
                elif xdb == "Reactome":
                    (xcols, xvals) = self._extract_cat_reactome(rconts, xpats)
                elif xdb == "CORUM":
                    (xcols, xvals) = self._extract_cat_corum(corum_json, xpats, acc)
                elif xdb == "MIM":
                    (xcols, xvals) = self._extract_cat_mim(record.comments, rconts, xpats)
                elif xdb == "DrugBank":
                    (xcols, xvals) = self._extract_cat_drugbank(rconts, xpats)
            else:
                # create empty dict with the name of colum
                for xc,xv in xpats:
                    xcols.append(xc)
                    xvals.append([np.nan])
                xvals = list(map(list, zip(*xvals)))
            # create dataframe with the Xreferece information
            df2 = pd.DataFrame(columns=xcols, data=xvals)
            if not df2.dropna(how='all').empty:
                # we add the given Isoforms ids
                dfx = pd.DataFrame(columns=[self.XID], data=IsoIds)
                df2 = pd.concat([df2,dfx],axis=1).ffill()
                df2.set_index(self.XID, inplace=True)
                # join using the index which is the UniProt accession of isoform
                df1 = df1.join(df2, how='outer')


        if 'xref_Ensembl_GeneId' in df1.columns and 'xref_Ensembl_transcId' in df1.columns:
            # reset index
            df1 = df1.reset_index()
            
            # extract the APPRIS decision and MANE from cross-references data ---
            if not appris_df.empty:
                df1 = df1.merge(appris_df, how='left', left_on=['xref_Ensembl_GeneId','xref_Ensembl_transcId'], right_on=['Gene ID','Transcript ID']).drop_duplicates()
                df1 = df1.drop(columns=['Gene ID','Transcript ID'])
    
            # extract the TRIFID score from cross-references data ---
            if not trifid_df.empty:
                df1 = df1.merge(trifid_df, how='left', left_on=['xref_Ensembl_GeneId','xref_Ensembl_transcId'], right_on=['gene_id','transcript_id']).drop_duplicates()
                df1 = df1.drop(columns=['gene_id','transcript_id'])
    
            # extract the CORSAIR score from cross-references data ---
            if not corsair_df.empty:
                df1 = df1.merge(corsair_df, how='left', left_on=['xref_Ensembl_GeneId','xref_Ensembl_transcId'], right_on=['ensembl_gene_id','ensembl_transc_id']).drop_duplicates()
                df1 = df1.drop(columns=['ensembl_gene_id','ensembl_transc_id'])
            
            # set index (Protein)
            df1.set_index(self.XID, inplace=True)

        return df1
    
    def _extract_cat_go(self, rconts, xpats):
        '''
        Parse the GO record
        '''
        xcols = []
        xvals = []
        # create a dict based on the GO type and the evidence codes as keys
        # rconts = [('GO:0016021', 'C:integral component of membrane', 'IEA:UniProtKB-KW')
        # filters the GO terms by:        
        flts = ['EXP','IDA','IPI','IMP','IGI','IEP','HTP','HDA','HMP','HGI','HEP','IBA','IBD','IKR','IRD']
        rcs = {}
        for rcont in rconts:
            g = rcont[1].split(':')[0]
            c = rcont[2].split(':')[0]
            if c in flts:
                # s = rcont[0]+'>'+rcont[1]+'|'+rcont[2]
                s = rcont[0]+'>'+rcont[1] # don't display the evidence code
                if g in rcs:
                    # rcs[g] += f";{s}"
                    rcs[g] += f"//{s}"
                else:
                    rcs[g] = s
        # go through all columns of xterms
        for xpat in xpats:
            xc = xpat[0] # column name
            g = xpat[1] # go type
            # create list of cols and values
            if g in rcs:
                r = rcs[g]
                xcols.append(xc)
                xvals.append([r])
            else:
                xcols.append(xc)
                xvals.append([np.nan])
        # transpose list of lists
        xvals = list(map(list, zip(*xvals)))
        return (xcols, xvals)

    def _extract_cat_kegg(self, datadict, rconts, xpats):
        '''
        Parse the KEGG record
        '''
        xcols = []
        xvals = []
        # go through all columns of xterms
        for xpat in xpats:
            xc = xpat[0] # column name
            # xp = xpat[1] # pattern
            # extract the record value that achives the pattern
            rcs = []
            for rcont in rconts:
                id = rcont[0]
                rc = ''
                try:
                    record = datadict[id]
                    if record:
                        pattern = re.search(r'(PATHWAY\s*[\w\W]*)', record, re.I | re.M)
                        if pattern:
                            for m in pattern[1].split('\n'):
                                if m.startswith('PATHWAY'):
                                    m = re.sub(r'PATHWAY\s*','',m).strip()
                                    ms = re.split(r'\s+', m, 1) # split only for the first space
                                    # rc += f"{ms[0]}>{''.join(ms[1:])};"
                                    rc += f"{ms[0]}>{''.join(ms[1:])}//"
                                elif m.startswith(' '):
                                    m = re.sub(r'^\s*','',m).strip()
                                    ms = re.split(r'\s+', m, 1) # split only for the first space
                                    # rc += f"{ms[0]}>{''.join(ms[1:])};"
                                    rc += f"{ms[0]}>{''.join(ms[1:])}//"
                                else:
                                    break
                        if rc != '':
                            # rc = re.sub(r'\;$','', rc) # delete ; at the end of string
                            rc = re.sub(r'\/\/$','', rc)
                    if rc != '': rcs.append(rc)
                        
                    pass
                except:
                    pass
            # create list of cols and values
            if rcs:
                # rcs = ";".join(rcs)
                rcs = "//".join(rcs)
                xcols.append(xc)
                xvals.append([rcs])
            else:
                xcols.append(xc)
                xvals.append([np.nan])
        # transpose list of lists
        xvals = list(map(list, zip(*xvals)))
        return (xcols, xvals)

    def _extract_cat_panther(self, df, rconts, xpats, acc):
        '''
        Parse the PANTHER record
        '''
        xcols = []
        xvals = []
        # go through all columns of xterms
        for xpat in xpats:
            xc = xpat[0] # column name            
            try:
                # extract the information from the Panther id
                rcs = ''
                for rcont in rconts:
                    id = rcont[0]
                    x = df[df[3].str.startswith(id)][[3,4]].values.tolist()[0] # get the panther id and family description
                    x[0] = re.sub(r'\:.*$','',x[0]) # remove the subfamily id
                    # dsc = x[1].replace(';',',')
                    # rcs = f"{x[0]}>{dsc};"
                    rcs = f"{x[0]}>{x[1]}//"
                pass
            except Exception:
                rcs = ''
                pass            
            # create list of cols and values
            if rcs != '':
                rcs = re.sub(r'\;$','', rcs) # delete ; at the end of string
                rcs = re.sub(r'\/\/$','', rcs)
                xcols.append(xc)
                xvals.append([rcs])
            else:
                xcols.append(xc)
                xvals.append([np.nan])
        # transpose list of lists
        xvals = list(map(list, zip(*xvals)))
        return (xcols, xvals)

    def _extract_cat_reactome(self, rconts, xpats):
        '''
        Parse the Reactome record
        '''
        xcols = []
        xvals = []
        # go through all columns of xterms
        for xpat in xpats:
            xc = xpat[0] # column name
            # xp = xpat[1] # pattern
            # extract the record value that achives the pattern
            rcs = []
            for rcont in rconts:
                id = rcont[0]
                dsc = "|".join(rcont[1:])
                dsc = re.sub(r'\s*\[[^\]]*\]\s*$','',dsc)
                dsc = re.sub(r'\s*\.\s*$','',dsc)
                # dsc = dsc.replace(';',',')
                rc = f"{id}>{dsc}"
                rcs.append(rc)
            # create list of cols and values
            if rcs:
                # rcs = ";".join(rcs)
                rcs = "//".join(rcs)
                xcols.append(xc)
                xvals.append([rcs])
            else:
                xcols.append(xc)
                xvals.append([np.nan])
        # transpose list of lists
        xvals = list(map(list, zip(*xvals)))
        return (xcols, xvals)

    def _extract_cat_corum(self, datatxt, xpats, acc):
        '''
        Parse the Corum record
        '''
        xcols = []
        xvals = []
        # go through all columns of xterms
        for xpat in xpats:
            xc = xpat[0] # column name
            # xp = xpat[1] # pattern
            # extract the information from the given UniProt accession
            rcs = []
            if datatxt:
                comps = list(filter(lambda person: acc in [i["swissprot"]["uniprot_id"] for i in person['subunits']], datatxt))
                if comps:
                    for comp in comps:
                        if 'complex_id' in comp and 'complex_name' in comp:
                            rcs += [f"compID_{comp['complex_id']}>{comp['complex_name']}"]
                if rcs:
                    rcs = "\\".join(rcs)
                else:
                    rcs = ''
            # create list of cols and values
            if rcs != '':
                xcols.append(xc)
                xvals.append([rcs])
            else:
                xcols.append(xc)
                xvals.append([np.nan])
        # transpose list of lists
        xvals = list(map(list, zip(*xvals)))
        return (xcols, xvals)

    def _extract_cat_mim(self, rcomms, rconts, xpats):
        '''
        Extract the MIM record with the description of the diseases
        '''
        xcols = []
        xvals = []
        # go through all columns of xterms
        for xpat in xpats:
            xc = xpat[0] # column name
            # xp = xpat[1] # pattern
            # extract the description of diseases from the UniProt comments
            rcs = ''
            if rcomms:
                # extract isoforms IDs and the displayed isoform
                rcomms = [c for c in rcomms if 'DISEASE:' in c]
                rcomms = [ re.findall(r"DISEASE:\s*([^\[]+)\[(MIM:\d+)\]\s*:", c) for c in rcomms ]
                rcomms = [ c[1]+'>'+c[0].strip() for rcom in rcomms for c in rcom if c ]
                # rcs += ";".join([ c.replace(';',',') for c in rcomms ])
                rcs += "//".join([ c for c in rcomms ])
            # create list of cols and values
            if rcs != '':
                xcols.append(xc)
                xvals.append([rcs])
            else:
                xcols.append(xc)
                xvals.append([np.nan])
        # transpose list of lists
        xvals = list(map(list, zip(*xvals)))
        return (xcols, xvals)

        
    def _extract_cat_drugbank(self, rconts, xpats):
        '''
        Parse the DrugBank record
        '''
        xcols = []
        xvals = []
        # go through all columns of xterms
        for xpat in xpats:
            xc = xpat[0] # column name
            # xp = xpat[1] # pattern
            # extract the record value that achives the pattern
            rcs = []
            for rcont in rconts:
                id = rcont[0]
                dsc = "|".join(rcont[1:])
                # dsc = dsc.replace(';',',')
                rc = f"{id}>{dsc}"
                rcs.append(rc)
            # create list of cols and values
            if rcs:
                # rcs = ";".join(rcs)
                rcs = "//".join(rcs)
                xcols.append(xc)
                xvals.append([rcs])
            else:
                xcols.append(xc)
                xvals.append([np.nan])
        # transpose list of lists
        xvals = list(map(list, zip(*xvals)))
        return (xcols, xvals)


    def to_file(self, df):
        '''
        Print to file
        '''
        df = df.reset_index()
        # add NaN values to the columns that do not exits in the current dataframe
        cols = df.columns.tolist()
        diff_cols = list(set(self.HEADERS) - set(cols))
        for c in diff_cols:
            df[c] = np.nan
        df.to_csv(self.outfile, sep="\t", index=False, columns=self.HEADERS)
        

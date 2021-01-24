import sys, os, logging
import urllib.request
import datetime
import re
import zipfile
import json
import pandas as pd
import numpy as np
from Bio import SwissProt
from Bio import SeqIO
from Bio.KEGG import REST


class creator:
    # https://www.uniprot.org/uniprot/?query=proteome:up000005640&format=fasta&include=yes&fil=reviewed:yes
    URL_UNIPROT = 'https://www.uniprot.org/uniprot/?'
    URL_UNIPROT += 'include=yes&' # include all isoforms
    URL_CORUM   = 'http://mips.helmholtz-muenchen.de/corum/download/allComplexes.json.zip'
    URL_PANTHER = 'ftp://ftp.pantherdb.org/sequence_classifications/current_release/PANTHER_Sequence_Classification_files/'
    SPECIES_LIST = {
        'human': {
            'scientific': 'Homo sapiens',
            'proteome': 'UP000005640'
        },
        'mouse': {
            'scientific': 'Mus musculus',
            'proteome': 'UP000000589'
        },
        'rat': {
            'scientific': 'Rattus norvegicus',
            'proteome': 'UP000002494'
        },
        'pig': {
            'scientific': 'Sus scrofa',
            'proteome': 'UP000008227'
        },
        'rabbit': {
            'scientific': 'Oryctolagus cuniculus',
            'proteome': 'UP000001811'
        },
        'zebrafish': {
            'scientific': 'Danio rerio',
            'proteome': 'UP000000437'
        }
    }
    # Column name with the cross-reference id
    XID = 'Isoform'
    # Meta terms of isoform
    META = ['xref_UniProt_Name','Isoform','Gene','prot_UniProt_Class','prot_Species','Description','Comment_Line']
    # Xreferences terms of isoform
    XTERMS = [
        ('Ensembl',  [('xref_Ensembl_protId','(ENS\w*P\d+[.]?\d*)'),('xref_Ensembl_transcId','(ENS\w*T\d+[.]?\d*)'),('xref_Ensembl_GeneId','(ENS\w*G\d+[.]?\d*)'),('Isoform','\[([^\]]*)\]')]),
        ('RefSeq',   [('xref_RefSeq_protId','(NP_\d+[.]?\d*)'),('xref_RefSeq_transcId','(NM_\d+[.]?\d*)'),('Isoform','\[([^\]]*)\]')]),
        ('CCDS',     [('xref_CCDS','(CCDS\d+[.]?\d*)'),('Isoform','\[([^\]]*)\]')]),
    ]
    # Category terms of isoform
    CTERMS = [
        ('GO',       [('cat_GO_C','^C:'),('cat_GO_F','^F:'),('cat_GO_P','^P:')]),
        ('KEGG',     [('cat_KEGG','')]),
        ('PANTHER',  [('cat_PANTHER','')]),
        ('Reactome', [('cat_Reactome','')]),
        ('CORUM',    [('cat_CORUM','')]),
        ('DrugBank', [('cat_DrugBank','')])
    ]
    HEADERS = [ h for h in META] + [ h[0] for i in XTERMS for h in i[1] if h[0] != 'Isoform' ] + [ h[0] for i in CTERMS for h in i[1] if h[0] != 'xref_UniProt_Acc' ]
    TIME = datetime.datetime.now().strftime("%Y%m")

    '''
    Creates the databases
    '''
    def __init__(self, s, o, f=None, d=False):
        
        # assign species
        species = s.lower()
        if species in self.SPECIES_LIST:
            self.species = species
            self.proteome_id = self.SPECIES_LIST[self.species]['proteome']
        else:
            sys.exit( "ERROR: Species parameter has been not found. Try with: "+", ".join(self.SPECIES_LIST.keys()) )
        
        # create output directory if does not exist
        self.outdir = o
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir, exist_ok=True)
            
        # create temporal file
        self.TMP_DIR = os.path.dirname(os.path.abspath(__file__)) +'/../tmp/'+ self.TIME +'/'+ self.species
        os.makedirs(self.TMP_DIR, exist_ok=True)
        logging.debug(self.TMP_DIR)
        
        # download sequences from UniProt
        self.outfname = species +'_'+ self.TIME +'_'+ f if f else ''
        self.db_fasta = self.outdir +'/'+ self.outfname +'.fasta'                
        self._download_fasta_db(self.db_fasta, f)
        
        # remove duplicate sequences
        if d:
            logging.debug("remove duplicate sequences")
            self._remove_duplicates(self.db_fasta, self.db_fasta)
            
        # create report from fasta file with the UniProt accession!! (key_function)
        self.db_fasta_seqio = SeqIO.index(self.db_fasta, "fasta", key_function=lambda rec : rec.split("|")[1])

        # create data files
        self.db_uniprot = self.TMP_DIR +'/'+ self.outfname +'.uniprot.dat'
        self.db_corum   = self.TMP_DIR +'/'+ ".".join(os.path.basename( self.URL_CORUM ).split(".")[:-1]) # get the filename from the URL (without 'zip' extension)
        self.db_panther = self.TMP_DIR +'/'+ self.outfname +'.panther.dat'
        
        # create output files
        self.outfile = self.outdir +'/'+ self.outfname +'.categories.tsv'


    def _delete_tmp_dir(self, dir):
        files = [ f for f in os.listdir(dir) ]
        for f in files:
            try:
                os.remove(os.path.join(dir, f))
            except Exception as e:
                logging.error(e)


    def _download_fasta_db(self, outfile, filt):
        '''
        Download the fasta database file
        '''
        url = self.URL_UNIPROT +'query=proteome:'+ self.proteome_id        
        if filt and filt == "sw": # filter by SwissProt
            url += '&fil=reviewed:yes'            
        elif filt and filt == "tr": # filter by TrEMBL
            url += '&fil=reviewed:no'
        url += '&format=fasta'
        logging.debug('get '+url)
        urllib.request.urlretrieve(url, outfile)


    def _remove_duplicates(self, infile, outfile=None):
        '''
        Remove duplicated sequences
        '''
        seqs = dict()
        precords = SeqIO.parse(infile, "fasta")
        records = SeqIO.to_dict(precords)
        # read all sequences
        # delete the duplicated sequences based on the sorted id's
        for i in list(records):
            record = records[i]
            s = str(record.seq)
            if s in seqs:
                l = [seqs[s], i]
                logging.warning( "duplicated sequences: {}".format(",".join(l)) , exc_info=False)
                if outfile:
                    l.sort()
                    logging.warning( "deleting the sequences {}".format(l[1]) , exc_info=False)
                    del records[l[1]]
                    seqs[s] = l[0]
            else:
                seqs[s] = i
        # write file if apply
        if outfile:
            with open(outfile, 'w') as handle:
                SeqIO.write(records.values(), handle, 'fasta')
        return None


    def download_raw_dbs(self, filt):
        '''
        Download the raw databases
        '''
        # delete any temporal file
        # self._delete_tmp_dir(self.TMP_DIR)
        
        # UniProt
        # filter by SwissProt (Reviewd) if apply
        if not os.path.isfile(self.db_uniprot):
            url = self.URL_UNIPROT +'query=proteome:'+ self.proteome_id
            if filt and filt == "sw":
                url += '&fil=reviewed:yes'
            url += '&format=txt'
            logging.debug("get "+url)
            urllib.request.urlretrieve(url, self.db_uniprot)
        else:
            logging.debug('cached uniprot')
        
        # CORUM
        # download all complexes file (using the same name)
        # unzip the file
        if not os.path.isfile(self.db_corum):
            url = self.URL_CORUM
            db_dat = self.TMP_DIR +'/'+ os.path.basename(url)
            logging.debug("get "+url)
            urllib.request.urlretrieve(url, db_dat)
            zip_ref = zipfile.ZipFile(db_dat, 'r')
            zip_ref.extractall(self.TMP_DIR)
            zip_ref.close()
        else:
            logging.debug('cached corum')
        
        # PANTHER
        # get the list of species and extract the file name
        if not os.path.isfile(self.db_panther):
            url = self.URL_PANTHER
            result = urllib.request.urlopen(url).read().decode('utf-8')
            if result:
                pattern = re.search(r'\s*(PTHR[^\_]*\_'+self.species+')', result, re.I | re.M)
                if pattern:
                    url = self.URL_PANTHER + pattern[1]
                    logging.debug("get "+url)
                    urllib.request.urlretrieve(url, self.db_panther)
        else:
            logging.debug('cached panther')
    
    
    def create_qreport(self):
        '''
        Create protein report
        '''
        # declare output dataframe
        df = pd.DataFrame()
        
        if self.db_uniprot:
            # create reports from external data ---
            logging.info('create reports from external data...')
            corum_json = None
            panther_txt = None
            if os.path.isfile(self.db_corum):
                with open(self.db_corum, 'r') as f:
                    corum_json = json.load(f)
            logging.debug('corum done')
            if os.path.isfile(self.db_panther):
                with open(self.db_panther, 'r') as f:
                    panther_txt = f.read()
            logging.debug('panther done')
            
            
            # Extract the info from the main database (UniProt), if apply
            # create cross-references data
            logging.info('create cross-references data from UniProtKB database...')
            for record in SwissProt.parse( open(self.db_uniprot) ):
                
                # extract main info ---
                name = record.entry_name
                acc = record.accessions[0]
                # accs = ";".join(record.accessions[1:])
                pattern = re.search(r'Name=(\w*)', record.gene_name, re.I | re.M)
                gene = pattern[1] if pattern else record.gene_name  
                pattern = re.search(r'[RecName|SubName]: Full=([^\;|\{]*)', record.description, re.I | re.M)
                dsc = pattern[1] if pattern else record.description
                dclass = record.data_class
                pattern = re.search(r'([\w|\s]*)\s+\(\w*\)', record.organism, re.I | re.M)
                species = pattern[1] if pattern else record.organism
                # extract isoforms IDs
                altprod = [c for c in record.comments if 'ALTERNATIVE PRODUCTS:' in c]
                if altprod:
                    IsoIds = re.findall(r'IsoId=([^\;|\,]*)', altprod[0], re.I | re.M | re.DOTALL)
                    # delete *-1 prefix from isoform Ids
                    IsoIds = [ i.replace('-1','') for i in IsoIds ]
                else:
                    IsoIds = [acc]
                # extract the comment line of fasta
                comm = ''
                if acc in self.db_fasta_seqio:
                    comm = ">"+self.db_fasta_seqio[acc].description
                    
                
                # create a dataframe with the Metadata information ---
                # UniProt accesion isoform is the index
                df1 = pd.DataFrame(columns=self.META, data=[[name,IsoIds,gene,dclass,species,dsc,comm]])
                df1 = df1.explode(self.XID)
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
                    xpats = terms[1]
                    xcols,xvals = [],[]
                    if xdb in rcross:
                        rconts = rcross[xdb]
                        if   xdb == "Ensembl":
                            (xcols, xvals) = self._extract_xref_ids(rconts, xpats, acc)
                        elif xdb == "RefSeq":
                            (xcols, xvals) = self._extract_xref_ids(rconts, xpats, acc)
                        elif xdb == "CCDS":
                            (xcols, xvals) = self._extract_xref_ids(rconts, xpats, acc)
                    else:
                        # create empty dict with the name of colum
                        for xc,xv in xpats:
                            xcols.append(xc)
                            xvals.append([np.nan])
                        xvals = list(map(list, zip(*xvals)))
                    # create dataframe with the Xreferece information
                    df2 = pd.DataFrame(columns=xcols, data=xvals)
                    if not df2.dropna().empty:
                        # check if UniProt accession does not exit
                        # we add the given Isoforms ids
                        if df2[self.XID].dropna().empty:
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

                # extract the category data ---
                for terms in self.CTERMS:
                    xdb = terms[0]
                    xpats = terms[1]
                    xcols,xvals = [],[]
                    if xdb in rcross:
                        rconts = rcross[xdb]
                        if xdb == "GO":
                            (xcols, xvals) = self._extract_cat_go(rconts, xpats)
                        elif xdb == "KEGG": # remote access
                            (xcols, xvals) = self._extract_cat_kegg(rconts, xpats)
                        elif xdb == "PANTHER":
                            (xcols, xvals) = self._extract_cat_panther(panther_txt, rconts, xpats, acc)
                        elif xdb == "Reactome":
                            (xcols, xvals) = self._extract_cat_reactome(rconts, xpats)
                        elif xdb == "CORUM":
                            (xcols, xvals) = self._extract_cat_corum(corum_json, xpats, acc)
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
                    if not df2.dropna().empty:
                        # we add the given Isoforms ids
                        dfx = pd.DataFrame(columns=[self.XID], data=IsoIds)
                        df2 = pd.concat([df2,dfx],axis=1).ffill()
                        df2.set_index(self.XID, inplace=True)
                        # join using the index which is the UniProt accession of isoform
                        df1 = df1.join(df2, how='outer')
                        
                
                # concatenate isoforms information
                df = pd.concat([df,df1])
                
                
        return df
    
    def _extract_xref_ids(self, rconts, xpats, acc):
        '''
        Parse the xref data
        '''
        xcols = []
        xvals = []
        # go through all columns of xterms
        for xpat in xpats:
            xc = xpat[0] # column name
            xp = xpat[1] # pattern
            # extract the record value that achives the pattern
            rcs = []
            for rcont in rconts:
                rc = [ re.findall(rf"{xp}", c) for c in rcont ]
                rc = "".join([i for s in rc for i in s]) # list of list to str
                # delete the last dot if apply
                rc = re.sub(r'\.$','', rc) if rc.endswith('.') else rc
                # exception
                rc = rc.replace('-1', '') if xc == self.XID else rc
                # only add values if exists
                if rc != '':
                    rcs.append(rc)
            # create list of cols and values
            if rcs:
                xcols.append(xc)
                xvals.append(rcs)
            else:
                xcols.append(xc)
                xvals.append([np.nan])
        # transpose list of lists
        xvals = list(map(list, zip(*xvals)))
        return (xcols, xvals)

    def _extract_cat_go(self, rconts, xpats):
        '''
        Parse the xref data
        '''
        xcols = []
        xvals = []
        # go through all columns of xterms
        for xpat in xpats:
            xc = xpat[0] # column name
            xp = xpat[1] # pattern
            # extract the record value that achives the pattern
            rcs = []
            for rcont in rconts:
                rc = [ rcont for c in rcont if re.search(rf"{xp}", c) ]
                if rc:
                    rc = rc[0] # first elem of comprehension list
                    rcs.append(rc)
            # create list of cols and values
            if rcs:
                rcs = ";".join([f"{c[0]}>{c[1]}|{c[2]}" for c in rcs])
                xcols.append(xc)
                xvals.append([rcs])
            else:
                xcols.append(xc)
                xvals.append([np.nan])
        # transpose list of lists
        xvals = list(map(list, zip(*xvals)))
        return (xcols, xvals)

    def _extract_cat_kegg(self, rconts, xpats):
        '''
        Parse the raw database file
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
                    record = REST.kegg_get(id).read()
                    if record:
                        pattern = re.search(r'DEFINITION\s*([^\n]*)', record, re.I | re.M)
                        rc += pattern[1] if pattern else ''
                        pattern = re.search(r'PATHWAY\s*([\w\W]*)MODULE', record, re.I | re.M)
                        rc += "|"+re.sub(r'\s*\n\s*','|', pattern[1]) if pattern else ''
                        rc = re.sub(r'[-|\|]*\s*$','', rc) # delete - or | at the end of string
                        rc = f"{id}>{rc}"
                        rcs.append(rc)
                    pass
                except:
                    pass
            # create list of cols and values
            if rcs:
                rcs = ";".join(rcs)
                xcols.append(xc)
                xvals.append([rcs])
            else:
                xcols.append(xc)
                xvals.append([np.nan])
        # transpose list of lists
        xvals = list(map(list, zip(*xvals)))
        return (xcols, xvals)

    def _extract_cat_panther(self, datatxt, rconts, xpats, acc):
        '''
        Parse the raw database file
        '''
        xcols = []
        xvals = []
        # go through all columns of xterms
        for xpat in xpats:
            xc = xpat[0] # column name
            # xp = xpat[1] # pattern
            # extract the information from the given UniProt accession
            rcs = ''
            if datatxt:
                pattern = re.search(rf"UniProtKB={acc}\t*([^\t]*)\t*([^\t]*)", datatxt, re.I | re.M)
                rcs += f"{pattern[1]}>{pattern[2]}" if pattern else ''
            # otherwise, we use the information from given records
            if rcs == '':
                rcs = ";".join([x[0] for x in rconts])
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

    def _extract_cat_reactome(self, rconts, xpats):
        '''
        Parse the rcont data
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
                rc = f"{id}>{dsc}"
                rcs.append(rc)
            # create list of cols and values
            if rcs:
                rcs = ";".join(rcs)
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
        Parse the raw database file
        '''
        xcols = []
        xvals = []
        # go through all columns of xterms
        for xpat in xpats:
            xc = xpat[0] # column name
            # xp = xpat[1] # pattern
            # extract the information from the given UniProt accession
            rcs = ''
            if datatxt:
                comps = list(filter(lambda person: acc in person['subunits(UniProt IDs)'], datatxt))
                if comps:
                    rcs += ";".join([ 'compID_'+str(comp['ComplexID'])+'>'+comp['ComplexName'] for comp in comps if 'ComplexID' in comp and 'ComplexName' in comp ])
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
        Parse the rcont data
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
                rc = f"{id}>{dsc}"
                rcs.append(rc)
            # create list of cols and values
            if rcs:
                rcs = ";".join(rcs)
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
        

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
            'taxonomy': '9606',
            'proteome': 'UP000005640'
        },
        'mouse': {
            'scientific': 'Mus musculus',
            'taxonomy': '10090',
            'proteome': 'UP000000589'
        },
        'rat': {
            'scientific': 'Rattus norvegicus',
            'taxonomy': '10116',
            'proteome': 'UP000002494'
        },
        'pig': {
            'scientific': 'Sus scrofa',
            'taxonomy': '9823',
            'proteome': 'UP000008227'
        },
        'rabbit': {
            'scientific': 'Oryctolagus cuniculus',
            'taxonomy': '9986',
            'proteome': 'UP000001811'
        },
        'zebrafish': {
            'scientific': 'Danio rerio',
            'taxonomy': '7955',
            'proteome': 'UP000000437'
        },
        'ecoli': {
            'scientific': 'Escherichia coli',
            'taxonomy': '562',
            'proteome': 'UP000000558'
        }
    }
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
        ('DrugBank', [('cat_DrugBank','')])
    ]
    HEADERS = [ h for h in META] + [ h for i in XTERMS for h in i[1][:-1] ] + [ h[0] for i in CTERMS for h in i[1] ]
    TIME = datetime.datetime.now().strftime("%Y%m")


    '''
    Creates the databases
    '''
    def __init__(self, s, o, f=None):
        
        # assign species
        species = s.lower()
        if species in self.SPECIES_LIST:
            self.species = species
            self.proteome_id = self.SPECIES_LIST[self.species]['proteome']
            self.scientific = self.SPECIES_LIST[self.species]['scientific']
            self.taxonomy = self.SPECIES_LIST[self.species]['taxonomy']
        else:
            sys.exit( "ERROR: Species parameter has been not found. Try with: "+", ".join(self.SPECIES_LIST.keys()) )
        
        # prepare temporal file
        self.TMP_DIR = os.path.dirname(os.path.abspath(__file__)) +'/../tmp/'+ self.TIME +'/'+ self.species
        os.makedirs(self.TMP_DIR, exist_ok=True)
        logging.debug(f"TMP_DIR: {self.TMP_DIR}")
        
        # prepare cached directory
        self.CHD_DIR = os.path.dirname(os.path.abspath(__file__)) +'/../cached/'+ self.species
        os.makedirs(self.CHD_DIR, exist_ok=True)
        logging.debug(f"CACHED_DIR: {self.CHD_DIR}")

        # prepare the output directory if does not exist
        self.outdir = os.path.dirname(o)
        if not os.path.exists(self.outdir):
            os.makedirs(self.outdir, exist_ok=True)
        
        # define output file
        self.outfile = o

        # define output files for "create_sb"
        self.db_uniprot = self.TMP_DIR +'/uniprot.dat'
        # self.db_uniprot = self.TMP_DIR +'/../../../test/test_1033.dat'
        
        self.db_fasta = self.TMP_DIR +'/proteins.fasta'
        self.db_corum   = self.TMP_DIR +'/'+ ".".join(os.path.basename( self.URL_CORUM ).split(".")[:-1]) # get the filename from the URL (without 'zip' extension)
        self.db_panther = self.TMP_DIR +'/panther.dat'
        self.cached_dir_kegg = self.CHD_DIR +'/kegg'
        os.makedirs(self.cached_dir_kegg, exist_ok=True)


    def _delete_tmp_dir(self, dir):
        files = [ f for f in os.listdir(dir) ]
        for f in files:
            try:
                os.remove(os.path.join(dir, f))
            except Exception as e:
                logging.error(e)

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


    def download_fasta_dbs(self, filt=None, d=None):
        '''
        Download the fasta file
        '''        
        # filter by SwissProt (Reviwed) if apply
        if not os.path.isfile(self.outfile):
            if filt and filt.startswith("pro"): # filter by proteome
                url = self.URL_UNIPROT +'query=proteome:'+ self.proteome_id
            else: # by default filter by organism
                url = self.URL_UNIPROT +'query='+ 'organism:'+self.species+ '+taxonomy:'+self.taxonomy
            if filt and filt.endswith("sw"): # filter by SwissProt
                url += '&fil=reviewed:yes'
            url += '&format=fasta'
            logging.info("get "+url+" > "+self.outfile)
            urllib.request.urlretrieve(url, self.outfile)
        else:
            logging.info('cached uniprot')
            
        # remove duplicate sequences
        if d:
            logging.info("remove duplicate sequences")
            self._remove_duplicates(self.outfile, self.outfile)
            
        
    def download_raw_dbs(self, filt=None):
        '''
        Download the raw databases
        '''        
        # UniProt Fasta
        if not os.path.isfile(self.db_fasta):
            url = self.URL_UNIPROT +'query='+ 'organism:'+self.species+ '+taxonomy:'+self.taxonomy
            url += '&format=fasta'
            logging.info("get "+url+" > "+self.db_fasta)
            urllib.request.urlretrieve(url, self.db_fasta)
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
            url = self.URL_UNIPROT +'query='+ 'organism:'+self.species+ '+taxonomy:'+self.taxonomy
            url += '&format=txt'
            logging.info("get "+url+" > "+self.db_uniprot)
            urllib.request.urlretrieve(url, self.db_uniprot)
        else:
            logging.info('cached uniprot data')
        
        # CORUM
        # download all complexes file (using the same name)
        # unzip the file
        if not os.path.isfile(self.db_corum):
            url = self.URL_CORUM
            db_dat = self.TMP_DIR +'/'+ os.path.basename(url)
            logging.info("get "+url+" > "+db_dat)
            urllib.request.urlretrieve(url, db_dat)
            zip_ref = zipfile.ZipFile(db_dat, 'r')
            zip_ref.extractall(self.TMP_DIR)
            zip_ref.close()
        else:
            logging.info('cached corum')
        
        # PANTHER
        # get the list of species and extract the file name
        if not os.path.isfile(self.db_panther):
            url = self.URL_PANTHER
            result = urllib.request.urlopen(url).read().decode('utf-8')
            if result:
                pattern = re.search(r'\s*(PTHR[^\_]*\_'+self.species+')', result, re.I | re.M)
                if pattern:
                    url = self.URL_PANTHER + pattern[1]
                    logging.info("get "+url+" > "+self.db_panther)
                    urllib.request.urlretrieve(url, self.db_panther)
        else:
            logging.info('cached panther')
    
    
    def create_qreport(self):
        '''
        Create protein report
        '''        
        if self.db_uniprot:
            # create reports from external data ---
            logging.info('create reports from external data...')
            corum_json = None
            panther_df = None
            if os.path.isfile(self.db_corum):
                with open(self.db_corum, 'r') as f:
                    corum_json = json.load(f)
            if os.path.isfile(self.db_panther):
                panther_df = pd.read_csv(self.db_panther, sep="\t", dtype=str, header=None, low_memory=False)
            
            logging.info("start the extraction of data...")
            ddf = []
            for rec in SwissProt.parse(open(self.db_uniprot)):
                ddf.append( self._create_qreport(rec,corum_json,panther_df) )
            ddf = pd.concat(ddf)

        return ddf


    def _create_qreport(self, record, corum_json, panther_df):
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
        pattern = re.search(r'Name=([^\s|\;]*)', record.gene_name, re.I | re.M)
        gene = pattern[1] if pattern else record.gene_name  
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
                        xvals.append([xp,xt,xg,i])                            
                elif xdb == "RefSeq":
                    for rcont in rconts:
                        xp = rcont[0]
                        # extract the isoform id from the last tuple (transc+isoform)
                        y = re.findall(r"\[([^\]]*)\]", rcont[1])
                        i = acc if not y or y[0] == IsoDisplayed else y[0]
                        # remove the isoform id from the last tuple
                        xt = re.sub(r"\.\s+.*$",'',rcont[1])                
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
                    (xcols, xvals) = self._extract_cat_panther(panther_df, rconts, xpats, acc)
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
                        
        return df1
    
    def _extract_cat_go(self, rconts, xpats):
        '''
        Parse the GO record
        '''
        xcols = []
        xvals = []
        # create a dict based on the GO type and the evidence codes as keys
        # rconts = [('GO:0016021', 'C:integral component of membrane', 'IEA:UniProtKB-KW')
        rcs = {}
        for rcont in rconts:
            g = rcont[1].split(':')[0]
            c = rcont[2].split(':')[0]
            if not g in rcs: rcs[g] = {}
            
            if c in rcs[g]:
                rcs[g][c].append(f"{rcont[0]}>{rcont[1]}|{rcont[2]}") 
            else:
                rcs[g][c] = [f"{rcont[0]}>{rcont[1]}|{rcont[2]}"]
        # go through all columns of xterms
        for xpat in xpats:
            xc = xpat[0] # column name
            g = xpat[1] # go type
            # create list of cols and values
            if g in rcs:
                z = rcs[g]
                r = ",".join([ f"{c}:[{';'.join(v)}]" for c,v in z.items() ])
                xcols.append(xc)
                xvals.append([r])
            else:
                xcols.append(xc)
                xvals.append([np.nan])
        # transpose list of lists
        xvals = list(map(list, zip(*xvals)))
        return (xcols, xvals)

    def _extract_cat_kegg(self, rconts, xpats):
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
                    id2 = id.replace(':','_')
                    of = self.cached_dir_kegg +f"/{id2}.txt"
                    of2 = self.cached_dir_kegg +f"/{id2}.dat"
                    if os.path.isfile(of):
                        with open(of, 'r') as f:
                            rc = f.read()
                    else:
                        if os.path.isfile(of2):
                            with open(of2, 'r') as f:
                                record = f.read()
                        else:
                            record = REST.kegg_get(id).read()
                            with open(of2, 'w') as f:
                                f.write(record)
                        if record:
                            pattern = re.search(r'DEFINITION\s*([^\n]*)', record, re.I | re.M)
                            rc += pattern[1] if pattern else ''
                            pattern = re.search(r'PATHWAY\s*([\w\W]*)MODULE', record, re.I | re.M)
                            rc += "|"+re.sub(r'\s*\n\s*','|', pattern[1]) if pattern else ''
                            rc = re.sub(r'[-|\|]*\s*$','', rc) # delete - or | at the end of string
                            rc = f"{id}>{rc}"
                            with open(of, 'w') as f:
                                f.write(rc)
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

    def _extract_cat_panther(self, df, rconts, xpats, acc):
        '''
        Parse the PANTHER record
        '''
        xcols = []
        xvals = []
        # go through all columns of xterms
        for xpat in xpats:
            xc = xpat[0] # column name
            # extract the information from the given UniProt accession
            try:
                x = df[df[1] == acc][[3,5]].values.tolist()[0]
                rcs = f"{x[0]}>{x[1]}"
                pass
            except Exception:
                rcs = ''
                pass            
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
        Parse the Corum record
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
        

#!/usr/bin/python

# import global modules
import os
import sys
import argparse
import logging
import re
import pandas as pd
import numpy as np


# Module metadata variables
__author__ = "Jose Rodriguez"
__credits__ = ["Jose Rodriguez", "Jesus Vazquez"]
__license__ = "Creative Commons Attribution-NonCommercial-NoDerivs 4.0 Unported License https://creativecommons.org/licenses/by-nc-nd/4.0/"
__version__ = "1.0.1"
__maintainer__ = "Jose Rodriguez"
__email__ = "jmrodriguezc@cnic.es"
__status__ = "Development"


###################
# Local functions #
###################

def get_cols_from_headers(cols, headers):
    out = []
    if cols and headers:
        if isinstance(headers, str): headers = re.split(r'\s*:\s*', headers)
        for c in headers:
            if c in cols:
                out.append(c)
            elif '*' in c:
                c = c.replace('*','\w+')
                for s in cols:
                    if re.match(c, s):
                        out = out + [ m for m in re.findall(c, s)]
            elif '[' in c and ']' in c:
                out.append(c)
            elif '{' in c and '}' in c:
                out.append(c)
    return out


def filter_rows(idf, filters):
    
    # filter values
    def _filter_tuple(tup, filt, tup_avail):
        # filter olny in the selected columns
        out = [ re.findall(rf"{filt}", str(t)) if i in tup_avail else t for i,t in enumerate(tup) ]
        # the list of tuples coming from the "findall" is joined with ';'
        # out = [ ';'.join([str(x) for t in l for x in t if x != '']) if isinstance(l,list) else l for l in out]
        out = [ '//'.join([str(x) for t in l for x in t if x != '']) if isinstance(l,list) else l for l in out]

        return tuple(out)
    
    # get columns
    cols = idf.columns.tolist()
    # get list of filters
    flts = re.split(r'\s*&\s*', filters) if filters else []
    for flt in flts:
        fc = re.split(r'\s*:\s*', flt)
        if fc and len(fc) == 2:
            col = fc[0]
            f = fc[1]            
            # get the columns that match with the given col name (regex)
            tup_avail = []
            col = col.replace('*','\w+')# convert to regex replacing '*' to '\w+'
            col = rf"({col})"
            for i,c in enumerate(cols):
                if re.match(col, c):
                    tup_avail = tup_avail + [ i for m in re.findall(col, c)]

            # if we find matched columns we apply the filters
            if tup_avail:
                # create pattern with the string separated by comma
                # eg.
                # f = IDA,ISS,HDA
                # r"IDA:\[([^\]]+)\]|ISS:\[([^\]]+)\]|HDA:\[([^\]]+)\]"
                p = ":\[([^\]]+)\]|".join( re.split(r'\s*,\s*',f) )
                p = rf"{p}:\[([^\]]+)\]"
                # create list of tuples
                # apply pattern for the list of tuples                    
                df_tuples = [tuple(x) for x in idf.to_numpy()]
                df_filt = [_filter_tuple(tup,p, tup_avail) for tup in df_tuples]
                idf = pd.DataFrame(df_filt, columns=idf.columns.tolist())
    
    # Important! for the reduction of memory
    # remove duplicates
    idf.drop_duplicates(inplace=True)

    return idf


def extract_and_join_columns(idf, headers_inf, headers_sup, cols_inf, cols_sup):
    # extract the columns. If there are multiple columns, join in one
    def _extract_and_join_columns(idf, cols, header):
        out = []
        if len(cols) > 1:
            if ':' in header:
                out = ["//".join([str(j) for j in s if not pd.isnull(j) and j != '']) for s in idf[cols].to_numpy()]
        elif len(cols) == 1:
            # get the column
            col = cols[0]
            # check if the col is a constant value: [X]
            # check if the col is a the order of column: {1}
            # otherwise, the col is the list of column names
            if col and '[' in col and ']' in col:
                c = re.findall(r'\[([^\]]*)\]', col)[0]
                # create a column with the given constant
                # extract the column with the constant
                idf[col] = c
                # idf = idf.reset_index()
                out = idf[col]
            elif col and '{' in col and '}' in col:
                c = int(re.findall(r'\{([^\}]*)\}', col)[0])
                # extract the column by position
                out = idf[idf.columns[c]]
            elif col:
                # extract the column
                out = idf[col]
        return out

    # create a list of tuple with the (input columns and the output heaers)
    colheaders = [(cols_inf,headers_inf)]
    if headers_sup and cols_sup: colheaders.append((cols_sup,headers_sup))
    
    # init output dataframe
    odf = pd.DataFrame(columns=[h[1] for h in colheaders])
    
    # get the columns with its -header
    for col,header in colheaders:
        odf[header] = _extract_and_join_columns(idf, col, header)

    # remove duplicates to reduce the use of memory
    odf.drop_duplicates(inplace=True)
    # remove row with any empty columns
    odf.replace('', np.nan, inplace=True)
    odf.dropna(inplace=True)
    
    return odf


def exploding_columns(idf):
    def _exploding_columns(idf, x, y):
        # replace np.nan to ''
        idf.replace(np.nan, '', inplace=True)
        # Exploding into multiple cells
        # We start with creating a new dataframe from the series with  as the index
        df = pd.DataFrame(idf[x].str.split('//').tolist(), index=idf[y].values).stack().rename(x)
        df = df.reset_index()
        # convert the index, which is a list of tuple, into columns
        # level_0 is the correct
        # remove level_1 column
        df.rename(columns={'level_0': y}, inplace=True)
        df.drop(['level_1'], axis=1, inplace=True)
        return df
        
    cols = idf.columns.tolist()
    df = idf
    for x in cols:
        # retrieve the column header other than the current loop
        y = "".join([i for i in cols if i != x])
        # check if '//' exits in column
        if any(df[x].str.contains('//')): df = _exploding_columns(df, x, y)
    return df


def replace_by_xrefprotein(intcols, iscols, df_inf, cols_datsup):
    '''
    WARNING!!! THIS METHOD IS HARD-CORE!! But it was the only way I thought to do it
    First, check if intersection columns is "Protein".
    Second, check if the columns of sup df contains the xref proteins.
    Third, extract the first value (not NaN) of 'Protein', and check if the value keeps the regex of one xref id.
    Fourth, replace the 'Protein' by the new xref column name.
    '''
    # check if intersection columns is "Protein"
    if 'Protein' in intcols:
        # check if the columns of sup df contains the xref proteins
        xrefs = {'xref_Ensembl_protId':   r'^ENS(\w{3})?P\d+',
                 'xref_Ensembl_transcId': r'^ENS(\w{3})?T\d+',
                 'xref_Ensembl_GeneId':   r'^ENS(\w{3})?G\d+',
                 'xref_RefSeq_protId':    r'^[N|X|Y]P_\d+',
                 'xref_RefSeq_transcId':  r'^[N|X][M|C]_\d+',
        }
        intxrefs = np.intersect1d(list(xrefs.keys()),cols_datsup).tolist() if cols_datsup else []
        # extract the first value (not NaN) of 'Protein'
        x = df_inf['Protein'].tolist()[0]
        # check if the value keeps the regex of one xref id
        m = [ i for i,r in xrefs.items() if re.match(r,x) ]
        if m:
            intcols = [ m[0] if i == 'Protein' else i for i in intcols ]
            iscols  = [ m[0] if i == 'Protein' else i for i in iscols ]
      
    return intcols,iscols

def merge_unknown_columns(df_inf, df_sup):
    # extract interseted column
    ic = df_inf.columns
    sc = df_sup.columns
    k = list(set(ic) & set(sc))
    if k and len(k) >= 1:
        df = pd.merge(df_inf, df_sup, on=k)
    else:
        df = None

    return df



#################
# Main function #
#################
def main(args):
    '''
    Main function
    '''
    # get input variables
    headers_inf = args.inf_headers
    headers_sup = args.sup_headers
    pattern = args.pattern

    
    logging.info("reading the input file...")
    outdat = pd.read_csv(args.infile, sep="\t", dtype=str, na_values=['NA'], low_memory=False)
    
    
    
    logging.info("getting the columns from the given headers...")
    all_cols = outdat.columns.to_list()
    cols_inf = get_cols_from_headers(all_cols, headers_inf)
    cols_sup = get_cols_from_headers(all_cols, headers_sup)
    
    
        

    logging.info("joining the columns...")
    outdat = extract_and_join_columns(outdat, headers_inf, headers_sup, cols_inf, cols_sup)
    
    

    logging.info("exploding the columns into multiple rows...")
    outdat = exploding_columns(outdat)


    logging.info("changing the order of columns...")
    cols = outdat.columns.to_list()
    cols = [cols[i] for i in [1,0,2] if (i < len(cols))]
    outdat = outdat[cols]

    
    # It is not necessary because the large file containing the categories no longer includes the evidence code.
    # pattern = r'\|[EXP|IDA|IPI|IMP|IGI|IEP|HTP|HDA|HMP|HGI|HEP|IBA|IBD|IKR|IRD].*$'
    if pattern:
        logging.info("remove a given pattern from the categories...")
        outdat[headers_inf].str.replace(pattern,'', regex=True, inplace=True)
        outdat[headers_sup].str.replace(pattern,'', regex=True, inplace=True)
        

    logging.info("remove duplicates and remove row with any empty column")
    # remove duplicates
    outdat.drop_duplicates(inplace=True)
    # remove row with any empty columns
    outdat.replace('', np.nan, inplace=True)
    outdat.dropna(inplace=True)

    logging.info('print output')
    outdat.to_csv(args.outfile, sep="\t", index=False)



if __name__ == "__main__":
    
    # parse arguments
    parser = argparse.ArgumentParser(
        description='Create a Relation Table (protein2category) based on the categories',
        epilog='''Examples:
* Create a Relation Table (protein2category) based on all the categories (cat_*)
    python src/create_rt.py    -ii databases/human_202206.categories.tsv -o databases/human_202206.q2c.tsv -i "Protein" -j "cat_*"

* Create a Relation Table (protein2category) based on the categories with the following headers: cat_GO_C, cat_GO_F, cat_GO_P, and cat_KEGG.
    python src/create_rt.py    -ii databases/human_202206.categories.tsv -o databases/human_202206.q2c.tsv -i "Protein" -j "cat_GO_C:cat_GO_F:cat_GO_P:cat_KEGG"
        ''',
        formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-ii', '--infile',  required=True, help='Input file')
    parser.add_argument('-i',  '--inf_headers',  required=True, help='Column(s) for the inferior level')
    parser.add_argument('-j',  '--sup_headers',  help='Column(s) for the superior level')
    parser.add_argument('-p',  '--pattern',  help='Regex pattern to remove from the category description')
    parser.add_argument('-o',  '--outfile', required=True, help='Output file with the relationship table')
    parser.add_argument('-vv', dest='verbose', action='store_true', help="Increase output verbosity")
    args = parser.parse_args()

    # get the name of script
    script_name = os.path.splitext( os.path.basename(__file__) )[0].upper()

    # logging debug level. By default, info level
    if args.verbose:
        logging.basicConfig(level=logging.DEBUG,
                            format=script_name+' - '+str(os.getpid())+' - %(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')
    else:
        logging.basicConfig(level=logging.INFO,
                            format=script_name+' - '+str(os.getpid())+' - %(asctime)s - %(levelname)s - %(message)s',
                            datefmt='%m/%d/%Y %I:%M:%S %p')

    # start main function
    logging.info('start script: '+"{0}".format(" ".join([x for x in sys.argv])))
    main(args)
    logging.info('end script')

# Scripts for the Databases in the Proteomics Unit of CNIC

### Requirements:
Python3 and the following Python packages:
- sys
- os
- logging
- urllib.request
- datetime
- re
- zipfile
- json
- from Bio import SwissProt
- from Bio import SeqIO
- from Bio.KEGG import REST


## Executions

### To execute the whole workflow:
```bash
./bin/create_db_sb.sh inprogress
```


### To execute the whole workflow given a specific version.
```bash
./bin/create_db_sb.sh 2
```
The above example will save the databases in the "{date}.2" folder


### The following scripts, download the FASTA sequences determining the source of database, UniProtKB or the UniProt Proteome:

#### Only SwissProt data
```bash
python src/create_fasta.py -s human  -f pro-sw         -o databases/human_202206_pro-sw.fasta     -vv  &> logs/create_fasta.human.log
```
#### With SwissProt+TrEMBL
```bash
python src/create_fasta.py -s human  -f pro-sw-tr      -o databases/human_202206_pro-sw-tr.fasta  -vv  &> logs/create_fasta.human.log
```
#### Remove duplicated sequences in the FASTA file based on the sorted id's
```bash
python src/create_fasta.py -s human  -f pro-sw-tr  -d  -o databases/human_202206_pro-sw-tr.fasta  -vv  &> logs/create_fasta.human.log
```


### The following scripts, create a system biology data with the whole UniProtKB for the given species:
```bash
python src/create_sb.py    -s human     -o databases/human_202206.categories.tsv       -vv  &> logs/create_sb.human.log
```


### The following scripts, create a relation table protein2category:

Create a Relation Table (protein2category) based on all the categories (cat_*)
```bash
python src/create_rt.py    -ii databases/human_202206.categories.tsv -o databases/human_202206.q2c.tsv -i "Protein" -j "cat_*"
```

Create a Relation Table (protein2category) based on the categories with the following headers: cat_GO_C, cat_GO_F, cat_GO_P, and cat_KEGG.
```bash
python src/create_rt.py    -ii databases/human_202206.categories.tsv -o databases/human_202206.q2c.tsv -i "Protein" -j "cat_GO_C:cat_GO_F:cat_GO_P:cat_KEGG"

```


# DecoyPYrat.v2 - Fast Hybrid Decoy Sequence Database Creation for Proteomic Mass Spectromtery Analyses

Improvement on DecoyPYrat script from Sanger (https://www.sanger.ac.uk/science/tools/decoypyrat)

## Execution

```bash
python src/decoyPYrat.v2.py  --output_fasta test/uniprot_MusMusculus_dic2016.decoy.fasta  --decoy_prefix=DECOY test/uniprot_MusMusculus_dic2016.fasta

cat uniprot_MusMusculus_dic2016.target.fa uniprot_MusMusculus_dic2016.decoy.fasta > uniprot_MusMusculus_dic2016.target-decoy.fa
```




# Scripts for the Databases in the Proteomics Unit of CNIC

### Requirements:
Python3 and the Python packages saved in the "requirements.txt" file


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

Create a Relation Table (protein2category) based on the categories with the following headers: cat_GO_C, cat_GO_F, cat_GO_P, cat_KEGG, cat_PANTHER and cat_Reactome.
```bash
python src/create_rt.py    -ii databases/human_202206.categories.tsv -o databases/human_202206.q2c.tsv -i "Protein" -j "cat_GO_C:cat_GO_F:cat_GO_P:cat_KEGG:cat_PANTHER:cat_Reactome"

```


# DecoyPYrat.v2 - Fast Hybrid Decoy Sequence Database Creation for Proteomic Mass Spectromtery Analyses

Improvement on DecoyPYrat script from Sanger (https://www.sanger.ac.uk/science/tools/decoypyrat)

## Execution

```bash
python src/decoyPYrat.v2.py  --output_fasta test/uniprot_MusMusculus_dic2016.decoy.fasta  --decoy_prefix=DECOY test/uniprot_MusMusculus_dic2016.fasta

cat uniprot_MusMusculus_dic2016.target.fa uniprot_MusMusculus_dic2016.decoy.fasta > uniprot_MusMusculus_dic2016.target-decoy.fa
```


# Add Human Orthologs

## Download the human orthologs from Ensembl Biomart

```bash
python src/download_orthologs.py --species 'rabbit' --outfile test/rabbit/biomart.tsv
```
The list of species are: (rat pig rabbit zebrafish chicken) ['rnorvegicus','sscrofa','ocuniculus','drerio','ggallus']



## Extract categories from orthologous genes


Retrieve the human categories from the orthologous genes:
```bash
python  src/categorize_orthologs.py  -im test/rabbit/biomart.tsv -ic1 test/human_202306.uniprot.tsv -ic2 test/rabbit/rabbit_202306.uniprot.tsv -o test/rabbit/rabbit_202306_with_human_orthologs.uniprot.tsv
```


Obtain the Relation Table with the human categories from the orthologous genes:

```
python src/create_rt.py   -ii test/rabbit/rabbit_202306_with_human_orthologs.uniprot.tsv -o test/rabbit/rabbit_202306_with_human_orthologs.q2c.tsv -i "Protein" -j "cat_GO_C:cat_GO_F:cat_GO_P:cat_KEGG:cat_PANTHER:cat_Reactome" -nj Categories
```

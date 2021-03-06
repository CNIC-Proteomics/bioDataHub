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

To execute the whole workflow:
```bash
./bin/create_db_sb.sh
```
To execute the whole workflow given a specific version.
```bash
./bin/create_db_sb.sh 2
```
The above example will save the databases in the "{date}.2" folder


The following scripts, download the FASTA sequences from UniProt proteomes and create a file with the system biology data for the given species:

Only SwissProt data
```bash
python src/create_db_sb.py -s human      -o databases -f sw -vv  &> logs/create_db_sb.human.log
python src/create_db_sb.py -s mouse      -o databases -f sw -vv  &> logs/create_db_sb.mouse.log
python src/create_db_sb.py -s pig        -o databases -f sw -vv  &> logs/create_db_sb.pig.log
python src/create_db_sb.py -s rabbit     -o databases -f sw -vv  &> logs/create_db_sb.rabbit.log
python src/create_db_sb.py -s zebrafish  -o databases -f sw -vv  &> logs/create_db_sb.zebrafish.log
```

With SwissProt+TrEMBL
```bash
python src/create_db_sb.py -s human      -o databases -vv  &> logs/create_db_sb.human.log
python src/create_db_sb.py -s mouse      -o databases -vv  &> logs/create_db_sb.mouse.log
python src/create_db_sb.py -s pig        -o databases -vv  &> logs/create_db_sb.pig.log
python src/create_db_sb.py -s rabbit     -o databases -vv  &> logs/create_db_sb.rabbit.log
python src/create_db_sb.py -s zebrafish  -o databases -vv  &> logs/create_db_sb.zebrafish.log
```

Remove duplicated sequences in the FASTA file based on the sorted id's
```bash
python src/create_db_sb.py -s human      -o databases -d -vv  &> logs/create_db_sb.human.log
python src/create_db_sb.py -s mouse      -o databases -d -vv  &> logs/create_db_sb.mouse.log
python src/create_db_sb.py -s pig        -o databases -d -vv  &> logs/create_db_sb.pig.log
python src/create_db_sb.py -s rabbit     -o databases -d -vv  &> logs/create_db_sb.rabbit.log
python src/create_db_sb.py -s zebrafish  -o databases -d -vv  &> logs/create_db_sb.zebrafish.log
```


Add into the crontab
```bash
crontab crontab/crontab_db_sb.sh
```

# DecoyPYrat.v2 - Fast Hybrid Decoy Sequence Database Creation for Proteomic Mass Spectromtery Analyses

Improvement on DecoyPYrat script from Sanger (https://www.sanger.ac.uk/science/tools/decoypyrat)

## Execution

```bash
python src/decoyPYrat.v2.py  --output_fasta test/uniprot_MusMusculus_dic2016.decoy.fasta  --decoy_prefix=DECOY test/uniprot_MusMusculus_dic2016.fasta

cat uniprot_MusMusculus_dic2016.target.fa uniprot_MusMusculus_dic2016.decoy.fasta > uniprot_MusMusculus_dic2016.target-decoy.fa
```




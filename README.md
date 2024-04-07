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



<br><br><br>
___

# DecoyPYrat.v2 - Fast Hybrid Decoy Sequence Database Creation for Proteomic Mass Spectromtery Analyses

Improvement on DecoyPYrat script from Sanger (https://www.sanger.ac.uk/science/tools/decoypyrat)

## Execution

```bash
python src/decoyPYrat.v2.py  --output_fasta test/uniprot_MusMusculus_dic2016.decoy.fasta  --decoy_prefix=DECOY test/uniprot_MusMusculus_dic2016.fasta

cat uniprot_MusMusculus_dic2016.target.fa uniprot_MusMusculus_dic2016.decoy.fasta > uniprot_MusMusculus_dic2016.target-decoy.fa
```



<br><br><br>
___

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



<br><br><br>
___

# Dowload AlphaFold proteomes


We have downloaded the AlphaFold PDB's from EBI:
https://alphafold.ebi.ac.uk/download#proteomes-section


To execute the whole workflow given a specific version.
```bash
./bin/download_af.sh -v v4 -f 202311
```
The above example will save the databases in the "{date}.v4" folder




<br><br><br>
___

# Calculate the Define Secondary Structure of Proteins (DSSP)


## Install DSSP

DSSP is the name of the original program, and since then multiple implementations have been developed. Following the dssp website , you can get the source code directly from github.
I am not sure how things work on windows, but you need to make sure that Biophyton points to the version that is installed in your system.

In ubuntu dssp is called mkdssp, so you need to
```
sudo apt-get install dssp
sudo ln -s /usr/bin/mkdssp /usr/bin/dssp
```

## Calculate DSSP

The DSSP program calculates the most likely secondary structure assignment by reading the position of the atoms in a protein followed by 
the calculation of the H-bond energy between all atoms.
https://biopython.org/docs/1.75/api/Bio.PDB.DSSP.html

The RSA of an amino acid residue is defined as the ratio of the solvent-exposed surface area of that residue observed in a given structure and the maximum obtainable value of the solvent-exposed surface area for this amino acid.
Thus, RSA adopts values between 0 and 1, with 0 corresponding to a fully buried and 1 to a fully accessible residue, respectively.
https://en.wikipedia.org/wiki/Relative_accessible_surface_area

You're correct in assuming that you need to choose a threshold below which a residue is considered "buried", or inaccessible.
However, you should not use the absolute accessible area; what you need to calculate is the relative solvent accessible area, or RSA.

To do this, you'll need to parse the output from DSSP, extract the value from the ACC column and divide it by the total surface area for 
the residue. You can find tables with surface areas of amino acids on the Web; here's one, or search Google.

Note that these areas are often estimated assuming that the residue, X, is in a tripeptide, G-X-G. This can result in estimates of RSA 
that are more than 100% (i.e. the area given by DSSP for highly-exposed residues can be more than the estimated value for the residue in G-X-G).

Once RSA is calculated, you need to select a threshold - and this is quite arbitrary. Some people define RSA < 20% as "completely buried".
It's probably best to read a few papers on the topic and see what most people use: this one is a good starting point (https://academic.oup.com/nar/article/33/10/3193/1009111).

Personally, I find the output of STRIDE somewhat easier to parse than DSSP (it does the same thing); also, you'll benefit from using an existing 
library for parsing (e.g. Bioperl's Bio::Structure::SecStr::DSSP::Res).


Calculate the DSSP printing the following info
The dssp data returned for a single residue is a tuple in the form:

|Tuple Index| Value             |
|:----------|:------------------|
| 0         |DSSP index         |
| 1         |Amino acid         |
| 2         |Secondary structure|
| 3         |Relative ASA       |
| 4         |Phi                |
| 5         |Psi                |
| 6         |NH–>O_1_relidx     |
| 7         |NH–>O_1_energy     |
| 8         |O–>NH_1_relidx     |
| 9         |O–>NH_1_energy     |
| 10        |NH–>O_2_relidx     |
| 11        |NH–>O_2_energy     |
| 12        |O–>NH_2_relidx     |
| 13        |O–>NH_2_energy     |


The DSSP codes for secondary structure used here are:

| Code | Structure                   |
|:-----|:----------------------------|
| H    | Alpha helix (4-12)          |
| B    | Isolated beta-bridge residue|
| E    | Strand                      |
| G    | 3-10 helix                  |
| I    | Pi helix                    |
| T    | Turn                        |
| S    | Bend                        |
| -    | None                        |


```bash
./bin/calc_dssp.sh -v v4 -w 50 -i /mnt/tierra/U_Proteomica/UNIDAD/Databases/AlphaFold -o /mnt/tierra/U_Proteomica/UNIDAD/Databases/DSSP
```




<br><br><br>
___

# Dowload DisProt regions

DisProt is the major manually curated repository of Intrinsically Disordered Proteins, both for structural and functional aspects. Expert curators are involved in collecting experimentally confirmed biological data, valuable for the scientific community, and for updating and maintaining DisProt over time. Learn more about DisProt and how to contribute to data curation biocuration page. Visit the browse section to explore the data collected in DisProt.

```bash
mkdir /mnt/tierra/U_Proteomica/UNIDAD/Databases/DisProt/2023_12
cd /mnt/tierra/U_Proteomica/UNIDAD/Databases/DisProt/2023_12
wget "https://disprot.org/api/search?release=2023_12&show_ambiguous=true&show_obsolete=false&format=tsv&namespace=all&get_consensus=false" -O DisProt_2023_12_ambiguous_evidences.tsv
```




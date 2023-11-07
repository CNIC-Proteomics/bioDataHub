___
## v2.12
```
DATE: 2023_11
```

+ Add the **Human Orthologs** for the species: 'rnorvegicus','cgriseus','sscrofa','ocuniculus','drerio','ggallus','btaurus'.

+ Add the label of Principal Isoform from APPRIS database, the TRIFID scores and CORSAIR score for the species: 'rat

+ The addition of the cRAP database has been disabled for the moment.

___
## v2.11
```
DATE: 2023_06
```

+ The file structure has been updated. Now, for each species, there are two folders: "sequences" and "categories".

+ The cRAP database (https://www.thegpm.org/crap/) has been added to the 'create_fasta' program.
Please note that the cRAP fasta has been modified to include the 'cRAP_' prefix for all proteins.

+ Now, we retrieve only the Proteome proteins.

+ The duplicated proteins from TrEMBL have been removed.

+ The 'Gallus gallus' (chicken) has been added.

+ The 'pdesc2cat' files have been discarded.

+ The 'create_rt' program creates the Relation Table protein2category by filtering the categories. This program is based on 'createRels.v0211.py'.

+ The large file containing the categories does not include the evidence code in GO.

+ The 'create_sb' program creates a report file with the following categories:
     GO, KEGG, PANTHER, Reactome, CORUM, MIM, and DrugBank.
However, the 'create_rt' program will filter the categories based on their type.

+ There is a single relation table for each category, as well as a relation table that includes the most common databases:
    GO, KEGG, PANTHER, and Reactome.


___
## v2.10
```
DATE: 2022_12
```

+ Add the species parameters from config file.

+ Catch the exception when PANTHER url does not exit.

___
## v2.9
```
DATE: 2022_11
```

+ Updating the URL's due the changes from PANTHER.

+ Updating the changes in the UniProt text report.

___
## v2.8
```
DATE: 2022_11
```

+ Updating the URL's due the changes from UniProtKB

___
## v2.7
```
DATE: 2022_06
```

+ The categories are separated by the "//" delimiter.

+ The program "createRels" separate the categories by "//".

+ Something was wrong in the decoyPyrat "fastas.close()" than after was correct.

+ Include a program that reports some statistic values.

___
## v2.6
```
DATE: 2021_10
```

+ And again come back to use the ";" as delimiter.

+ Include the program "createRels" in this repository for ever. It retrieves the Relation Tables proteinId2category and proteindescription2category.

___
## v2.5
```
DATE: 2021_08
```

+ The delimiter that divides the categories of a protein, has changed from ";" to "||".
___
## v2.4
```
DATE: 2021_05
```

+ Fixing some bugs in the parser section.

+ OMIM database has been included but then it has been disabled.

+ The following filter in the GO terms is by default:
"cat_GO_*:EXP,IDA,IPI,IMP,IGI,IEP,HTP,HDA,HMP,HGI,HEP,IBA,IBD,IKR,IRD"

+ The descriptions of PANTHER have increased.

___
## v2.3
```
DATE: 2021_05
```

+ DrugBank is disabled.

+ Retrieve the PANTHER id (without subfamily id) and the Family description

+ Retrieve the KEGG Pathways.

+ Include createRels programa version 0.2.11

___
## v2.2
```
DATE: 2021_05
```

+ Change the format GO categories based on the evidences codes.

+ Rewrite the obsolete code in the calculation of cross-references.

+ Remote the version in the identifiers (RefSeq and GENCODE)

___
## v2.1
```
DATE: 2021_04.1
```

+ Fixing a bug getting the gene name.

+ Fixing a bug getting the cross-references.

___
## v2.0
```
DATE: 2021_04
```

+ The program has been divided in two: create_fasta and create_sb
    - create_fasta, downloads the FASTA file.
    - create_sb, create the System Biology table.
    
+ Now you can download the proteins from the Proteome (pro) or the proteins from UniProtKB (uni).

+ The system downloads systematically the following databases: ("pro-sw" "pro-sw-tr" "uni-sw" "uni-sw-tr").

+ The protein length has been included in the record table.

+ The time of execution has been decreased!!

+ Fixing a problem with the displayed sequences and the cross-reference identifies.

+ The taxonomy and the name of species are the main parameter in the query of URL.

___
## v1.4
```
DATE: 2021_03
```

+ The output folder has changes. It includes a prefix per data. For example, 202103 (first), 202103.1, 202103.2, etc

+ The output folder has changed (//tierra.cnic.es/SC/U_Proteomica/UNIDAD/iSanXoT_DBs)

___
## v1.3
```
DATE: 2021_03
```

+ Bug fixed: the columns of categories do not correspond with the data.

___
## v1.2
```
DATE: 2021_02
```

+ We have included the default filter of "protein2category" file, 'cat_GO_*:EXP,IDA,IPI,IMP,IGI,IEP,HTP,HDA,HMP,HGI,HEP'

+ We have rename the columns: prot_Species -> Species and Isoform -> Protein

___
## v1.1
```
DATE: 2021_01
```

### Highlights

+ Now, it is reported the updated categories with old format of relationship file "protein2category".
+ The protein identifier has been included in the DECOY identifier of FASTA file.

### Changes in the Code

+ Add the "Comment_Line" columns in the "categories" file. With this columns, we can create the old relationship file "protein2category".
+ We have add the "createRels" program (from iSanXoT-0.2.3). This programs allow us to create the "protein2category" file.

___
## v1.0
```
DATE: 2020_10
```

### Highlights

+ The category file contains multiple columns with the UniProt information.
+ The date of the first files of this version is 202010.

### Changes in the Code

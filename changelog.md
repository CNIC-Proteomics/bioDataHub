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

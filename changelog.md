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

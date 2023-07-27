# Introduction

The programs for the Databases in the Proteomics Unit of CNIC are in the following repository:

https://github.com/CNIC-Proteomics/iSanXoT-dbscripts


In this folder are located the **protein sequences** and the **category files** for the proteins in multiple species:

- Human
- Mouse
- Rat
- Pig
- Rabbit
- Zebrafish
- E.Coli
- Chicken

# Description files

For each species, there are two folders: sequences and categories


## Sequences

This folder contains the protein sequences in FASTA format, which depend on the type of database.

The filename composition is:

    {species}_{date}_{database}.fasta

The type of *database* are the following:
- pro-sw    =>  SwissProt for the UniProtKB Proteome
- pro-sw-tr =>  SwissProt for the UniProtKB Proteome


There are target/decoy files for each type of protein sequence files described above.

The decoy files have been created by the DecoyPYrat script from Sanger (https://www.sanger.ac.uk/science/tools/decoypyrat)

The filename composition is:

    {species}_{date}_{database}.{target/decoy}.fasta

The type of *target/decoy* files are the following:
- target       =>  File with the protein sequences. In this file the peptides have been maked isobaric (replacing the I to L ).
- decoy        =>  File with the only the decoy sequences. In this file the decoys have been maked isobaric (replacing the I to L ).
- target-decoy =>  File that concatenate the target/decoy sequences (maked isobaric, replacing the I to L ).

## Categories

This folder contains the files containing functional information that categorizes the proteins.

The file structure name is:

    [q2c_|g2c_]{species}_{date}.{type_cat}.tsv

The files with the 'type_cat' suffix are as follows:

- *uniprot* "{species}_{date}.uniprot.tsv" => These files contain multiple columns with the following information:
    + Meta terms of Protein:
      * xref_UniProt_Name: UniProt Name
      * Protein
      * Gene: HUGO gene name
      * Species
      * Length: protein length
      * Description: UniProt protein description
      * Comment_Line: UniProt description that is within the FASTA file
      * prot_UniProt_Class:  If the entry belongs to the Swiss-Prot section of UniProtKB (reviewed) or to the computer-annotated TrEMBL section (unreviewed).
    + Xreferences terms of Protein
      * xref_Ensembl_protId: Ensembl Protein Id
      * xref_Ensembl_transcId: Ensembl Transcript Id
      * xref_Ensembl_GeneId: Ensembl Gene Id
      * xref_RefSeq_protId: RefSeq Protein Id
      * xref_RefSeq_transcId: RefSeq Trascript Id
      * xref_CCDS: CCDS Id
    + Category terms of Protein:
      * cat_GO_C: GO terms for Cellular component
      * cat_GO_F: GO terms for Molecular function
      * cat_GO_P: GO terms for Biological process
      * cat_KEGG: KEGG information
      * cat_PANTHER: PANTHER information
      * cat_Reactome: Reactome information
      * cat_CORUM: CORUM information
      * cat_OMIM: The OMIM database provides information on the disease(s) associated with genetic variations in a given protein
      * cat_DrugBank: Drug and drug target database

- *database name(s) for categories*: "q2c_{species}_{date}.{type_cat}.tsv"
This file represents a protein-to-category relation table, with the following columns:
    + The first column contains the functional description of a protein, as provided by the database named in the filename.
    + The second column contains the UniProt accession of the protein in question.

- *database name(s) for categories*: "g2c_{species}_{date}.{type_cat}.tsv"
This file represents a gene-to-category relation table, with the following columns:
    + The first column contains the functional description of a gene, as provided by the database named in the filename.
    + The second column contains the gene in question.


There is a statistical file, with the suffix 'stats', that displays the number of isoforms per category.


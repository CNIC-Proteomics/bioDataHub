#!/usr/bin/bash

# Constants
CODEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd -P)/.."
BASEDIR="/mnt/tierra/U_Proteomica/UNIDAD/Databases/APPRIS"
SPECIES_LIST=(human mouse rat pig zebrafish chicken)
# define UniProt database folder for the cross-reference identifiers
XREFDATE="202501"
XREFDIR="/mnt/tierra/U_Proteomica/UNIDAD/Databases/UniProt/${XREFDATE}"


# Control the parameters
if [[ ! -z "$1" ]]; then
  VERSION=".${1}"
else
  VERSION=''
fi


# Declare variables
DATE="$(date +"%Y%m")" # create date
OUTDIR="${BASEDIR}/${DATE}${VERSION}" # with date+version folder


# Function that executes the input command
run_cmd () {
  echo "-- $1"
  echo ""
  eval $1
}
# Concatenate species annotation files
process_file() {
  local infile="$1"
  local outfile="$2"
  local cmd=""
  if [ -f "${outfile}" ]; then
    cmd="tail -n +2 '${infile}' >> '${outfile}'"
  else
    cmd="cat '${infile}' > '${outfile}'"
  fi
  echo "${cmd}"
}


echo "going through the species..."
for SPECIES in "${SPECIES_LIST[@]}"
do
    # get local variables
    XREFFILE="${XREFDIR}/${SPECIES}/categories/${SPECIES}_${XREFDATE}.uniprot.tsv"
    METHODDIR="${OUTDIR}/${SPECIES}/tmp" # re-declare the outdir with date+version+species
    OUTFILE="${METHODDIR}/${SPECIES}_${DATE}.appris.tsv"

    # execute the program:
    # The following script downloads the annotations for the APPRIS methods that locate the annotation in a specific region of the protein.
    CMD1="python '${CODEDIR}/src/download_appris.py' -s ${SPECIES} -o '${METHODDIR}' -vv "
    # Convert the method annotations in GTF format to another GTF that references the protein region
    CMD2="python '${CODEDIR}/src/convert_appris.py' -ia '${METHODDIR}/*.gtf' -iu '${XREFFILE}' -o ${OUTFILE} -vv "
    run_cmd "${CMD1} && ${CMD2}"

done



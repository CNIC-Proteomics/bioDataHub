#!/usr/bin/bash

# Constants
CODEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd -P)/.."
BASEDIR="/mnt/tierra/U_Proteomica/UNIDAD/Databases/APPRIS"
BASESCRIPT="$(basename "$0")"
SPECIES_LIST=(human mouse rat pig zebrafish chicken)
METHOD="spade"


# Function to print usage information
print_usage() {
    echo "Usage: ${BASESCRIPT} -f <date> [options]"
    echo "Options:"
    echo "  -f  Date of annotations (ie. 202501) [REQUIRED]"
    echo "  -v  Version of annotations (ie. [4, inprogress])"
    echo "  -h  Display this help message"
}

# Control the parameters
while getopts "hf:v:w:" opt; do
    case $opt in
        f)
            DATE="${OPTARG}"
            ;;
        v)
            VERSION=".${OPTARG}"
            ;;
        h)
            print_usage
            exit 0
            ;;
        \?)
            echo "Invalid option: -${OPTARG}" >&2
            print_usage
            exit 1
            ;;
    esac
done
# check for required options
if [[ -z "${DATE}" ]]; then
    echo "Error: -f <date> option is required." >&2
    print_usage
    exit 1
fi



# Declare variables
# input folder/files
DATDIR="${BASEDIR}/${DATE}${VERSION}"
# output files
METHOD_CDS_GTF="${DATDIR}/spade.cds.gtf"
METHOD_PEP_GTF="${DATDIR}/spade.pep.gtf"
METHOD_TSV="${DATDIR}/spade.tsv"



# Function that executes the input command
run_cmd () {
  echo "-- $1"
  echo ""
  eval $1
}
# Concatenate files (without headers)
concat_files() {
  local infile="$1"
  local outfile="$2"
  local cmd="cat '${infile}' >> '${outfile}'"
  echo "${cmd}"
}
# Concatenate files (with header)
concat_files_wh() {
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


echo "## going through the species..."
for SPECIES in "${SPECIES_LIST[@]}"
do
    # get local variables
    INDIR_spe="${DATDIR}/${SPECIES}/tmp"
    OUTDIR_spe="${DATDIR}/${SPECIES}"
    FILE_CDS_GTF="${OUTDIR_spe}/${SPECIES}_${DATE}.spade.cds.gtf"
    FILE_PEP_GTF="${OUTDIR_spe}/${SPECIES}_${DATE}.spade.pep.gtf"
    FILE_TSV="${OUTDIR_spe}/${SPECIES}_${DATE}.spade.tsv"

    # create annotation file with the CDS coordinates in several formats
    CMD1="python '${CODEDIR}/src/add_cds_coords_spade.py' \
            -ia '${INDIR_spe}/appris_method.spade.gtf' \
            -oc ${FILE_CDS_GTF} \
            -op ${FILE_PEP_GTF} \
            -o  ${FILE_TSV} \
          -vv "
    run_cmd "${CMD1}"

    # concatenate files
    run_cmd "$(concat_files "${FILE_CDS_GTF}" "${METHOD_CDS_GTF}")"
    run_cmd "$(concat_files "${FILE_PEP_GTF}" "${METHOD_PEP_GTF}")"
    run_cmd "$(concat_files_wh "${FILE_TSV}" "${METHOD_TSV}")"

done



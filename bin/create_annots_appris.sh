#!/usr/bin/bash

# Constants
CODEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd -P)/.."
BASEDIR="/mnt/tierra/U_Proteomica/UNIDAD/Databases/APPRIS"
BASESCRIPT="$(basename "$0")"
SPECIES_LIST=(human mouse rat pig zebrafish chicken)
METHOD="appris"


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
METHOD_CDS_GTF="${DATDIR}/appris.cds.gtf"
METHOD_PEP_GTF="${DATDIR}/appris.pep.gtf"
METHOD_TSV="${DATDIR}/appris.tsv"



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
# Separate the Principal Isoform (PI) from the rest (Non-Principal Isoforms, NPI)
separate_PI() {
  local infile="$1"
  local pi_outfile="$2"
  local npi_outfile="$3"
  local cmd=""
  # check if the input file exists
  if [ -f "${infile}" ]; then
    # extract Principal Isoforms
    cmd1="grep PRINCIPAL '${infile}' > '${pi_outfile}'"
    # extract Non-Principal Isoforms
    cmd2="grep -v PRINCIPAL '${infile}' > '${npi_outfile}'"
    cmd="${cmd1} && ${cmd2}"
  else
    echo "Error: Input file ${infile} does not exist!" >&2
    return 1
  fi
  echo "${cmd}"
}
# Obtain the unique PI and NPI CDS coordiantes
intersect_1() {
  local infile_a="$1"
  local infile_b="$2"
  local outfile="$3"
  local cmd=""
  # check if the input file exists
  if [[ -f "${infile_a}" && -f "${infile_b}" ]]; then
    cmd="bedtools intersect -nonamecheck -a '${infile_a}' -b '${infile_b}' -wa -f 1 -v  >  '${outfile}'"
  else
    echo "Error: Input files: ['${infile_a}','${infile_b}'] do not exist!" >&2
    return 1
  fi
  echo "${cmd}"
}
# Obtain the intersection between coordinates using bedtools
intersect_2() {
  local infile_a="$1"
  local infile_b="$2"
  local outfile="$3"
  local cmd=""
  # check if the input file exists
  if [[ -f "${infile_a}" && -f "${infile_b}" ]]; then
    CMD1="bedtools intersect -nonamecheck -a '${infile_a}' -b '${infile_b}' -wa -f 1  >  '${outfile}'"
    CMD2="bedtools intersect -nonamecheck -a '${infile_b}' -b '${infile_a}' -wa -f 1  >> '${outfile}'"
    cmd="${CMD1} && ${CMD2}"
  else
    echo "Error: Input files: ['${infile_a}','${infile_b}'] do not exist!" >&2
    return 1
  fi
  echo "${cmd}"
}



echo "## going through the species..."
for SPECIES in "${SPECIES_LIST[@]}"
do
    # get local variables
    INDIR_spe="${DATDIR}/${SPECIES}/tmp"
    OUTDIR_spe="${DATDIR}/${SPECIES}"
    FILE_CDS_GTF="${OUTDIR_spe}/${SPECIES}_${DATE}.appris.cds.gtf"
    FILE_CDS_PI_GTF="${OUTDIR_spe}/${SPECIES}_${DATE}.appris.cds_pi.gtf"
    FILE_CDS_NPI_GTF="${OUTDIR_spe}/${SPECIES}_${DATE}.appris.cds_npi.gtf"
    FILE_CDS_INT_GTF="${OUTDIR_spe}/${SPECIES}_${DATE}.appris.cds_int.gtf"
    FILE_CDS_INT_PI_GTF="${OUTDIR_spe}/${SPECIES}_${DATE}.appris.cds_u_pi.gtf"
    FILE_CDS_INT_NPI_GTF="${OUTDIR_spe}/${SPECIES}_${DATE}.appris.cds_u_npi.gtf"
    FILE_CDS_OVR="${OUTDIR_spe}/${SPECIES}_${DATE}.appris.cds_overlap.gtf"
    FILE_PEP_GTF="${OUTDIR_spe}/${SPECIES}_${DATE}.appris.pep.gtf"
    FILE_TSV="${OUTDIR_spe}/${SPECIES}_${DATE}.appris.tsv"

    # create annotation file with the CDS coordinates in several formats
    CMD1="python '${CODEDIR}/src/add_cds_coords_appris.py' \
            -ia '${INDIR_spe}/appris_method.appris.gtf' \
            -ip '${INDIR_spe}/appris_data.pannot.gtf' \
            -oc ${FILE_CDS_GTF} \
          -vv "
    run_cmd "${CMD1}"

    # separate the Principal Isoform (PI) from the rest (Non-Principal Isoforms, NPI)
    run_cmd "$(separate_PI "${FILE_CDS_GTF}" "${FILE_CDS_PI_GTF}" "${FILE_CDS_NPI_GTF}")"

    # intersect the CDS coordinates
    run_cmd "$(intersect_1 "${FILE_CDS_PI_GTF}" "${FILE_CDS_NPI_GTF}" "${FILE_CDS_INT_PI_GTF}")"
    run_cmd "$(intersect_1 "${FILE_CDS_NPI_GTF}" "${FILE_CDS_PI_GTF}" "${FILE_CDS_INT_NPI_GTF}")"
    run_cmd "$(intersect_2 "${FILE_CDS_PI_GTF}" "${FILE_CDS_NPI_GTF}" "${FILE_CDS_INT_GTF}")"

    # create multiple files with the CDS coordinates in GTF format
    CMD1="python '${CODEDIR}/src/add_cds_intersect_appris.py' \
            -i  '${FILE_CDS_GTF}' \
            -ii '${FILE_CDS_INT_GTF}' \
            -iip '${FILE_CDS_INT_PI_GTF}' \
            -iin '${FILE_CDS_INT_NPI_GTF}' \
            -oc  '${FILE_CDS_OVR}' \
            -op  '${FILE_PEP_GTF}' \
            -ot  '${FILE_TSV}' \
          "
    run_cmd "${CMD1}"

    # concatenate files
    run_cmd "$(concat_files "${FILE_CDS_GTF}" "${METHOD_CDS_GTF}")"
    run_cmd "$(concat_files "${FILE_PEP_GTF}" "${METHOD_PEP_GTF}")"
    run_cmd "$(concat_files_wh "${FILE_TSV}" "${METHOD_TSV}")"

done



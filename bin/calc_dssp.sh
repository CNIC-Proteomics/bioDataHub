#!/usr/bin/bash

# Constants
CODEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd -P)/.."
BASESCRIPT="$(basename "$0")"
BASENAME="${BASESCRIPT%.*}"

    
# Function to print usage information
print_usage() {
    echo "Usage: ${BASESCRIPT} [options]"
    echo "Options:"
    echo "  -v  Version of AlphaFold (ie. v4)"
    echo "  -i  Input folder with AlphaFold predictions in PDB format"
    echo "  -o  Output folder with DSSP results"
    echo "  -w  Number of threads/n_workers"
    echo "  -h  Display this help message"
    # Add more options here if needed
}

# Control the parameters
while getopts "hv:i:o:w:" opt; do
    case $opt in
        v)
            VERSION="${OPTARG}"
            ;;
        i)
            INDIR="${OPTARG}"
            ;;
        o)
            OUTDIR="${OPTARG}"
            ;;
        w)
            NTHREADS="${OPTARG}"
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


# Declare variables
DATE="$(date +"%Y%m")" # create date
INDIR="${INDIR}/${DATE}.${VERSION}" # with date+version folder
OUTDIR="${OUTDIR}/${DATE}.${VERSION}" # with date+version folder
LOGDIR="${CODEDIR}/logs/${DATE}.${VERSION}" # with date+version folder
SPECIES_LIST=(human mouse rat zebrafish)

# Function that executes the input command
run_cmd () {
  echo "-- $1"
  echo ""
  eval $1
}

# prepare workspaces
mkdir "${LOGDIR}"



# CREATE FASTA files --------

# go through the species...
for SPECIES in "${SPECIES_LIST[@]}"
do
  # get local variables
  INDIR_spe="${INDIR}/${SPECIES}" # re-declare the outdir with date+version+species
  OUTDIR_spe="${OUTDIR}/${SPECIES}" # re-declare the outdir with date+version+species
  OUTNAME="${SPECIES}_${DATE}"
  LOGFILE="${LOGDIR}/${BASENAME}.${OUTNAME}.log"

  # execute commands
  CMD1="python '${CODEDIR}/src/calc_DSSP.py' -w ${NTHREADS} -i '${INDIR_spe}' -o '${OUTDIR_spe}' &> '${LOGFILE}' "
  run_cmd "${CMD1}"
done

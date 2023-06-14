#!/usr/bin/bash

# Declare variables
if [[ ! -z "$1" ]]; then
  VERSION=".${1}"
else
  VERSION=''
fi
DATE="$(date +"%Y%m")" # create date
# CODEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd -P)/.."
CODEDIR="S:/U_Proteomica/UNIDAD/DatosCrudos/jmrodriguezc/projects/iSanXoT-dbscripts"
BASEDIR="//tierra.cnic.es/SC/U_Proteomica/UNIDAD/iSanXoT_DBs"
OUTDIR="${BASEDIR}/${DATE}${VERSION}" # with date+version folder
WSDIR="${BASEDIR}/current_release"
LOGDIR="${CODEDIR}/logs/${DATE}${VERSION}" # with date+version folder

TYPE_LIST=("pro-sw" "pro-sw-tr")
SPECIES_LIST=(human mouse rat pig rabbit zebrafish ecoli chicken)

# Function that executes the input command
run_cmd () {
  echo "-- $1"
  echo ""
  eval $1
}

# prepare workspaces
mkdir "${LOGDIR}"

# CREATE FASTA files --------

# for the following databases and species...
for TYPE in "${TYPE_LIST[@]}"
do
  for SPECIES in "${SPECIES_LIST[@]}"
  do
    # get local variables
    OUTNAME="${SPECIES}_${DATE}_${TYPE}"
    OUTFILE="${OUTDIR}/${OUTNAME}.fasta"
    LOGFILE="${LOGDIR}/create_fasta.${OUTNAME}.log"

    OUTFILE_dc="${OUTDIR}/${OUTNAME}.decoy.fasta"
    OUTFILE_tg="${OUTDIR}/${OUTNAME}.target.fasta"
    OUTFILE_dc_tg="${OUTDIR}/${OUTNAME}.target-decoy.fasta"
    LOGFILE_dc_tg="${LOGDIR}/decoyPYrat.${OUTNAME}.log"

    # execute commands
    CMD1="python '${CODEDIR}/src/create_fasta.py' -s ${SPECIES} -f ${TYPE} -o '${OUTFILE}' -d -vv  &> '${LOGFILE}' "
    CMD2="python '${CODEDIR}/src/decoyPYrat.v2.py' --output_fasta '${OUTFILE_dc}' --decoy_prefix=DECOY -t '${OUTFILE}.tmp' '${OUTFILE}' &> '${LOGFILE_dc_tg}' && cat ${OUTFILE_tg} ${OUTFILE_dc} > ${OUTFILE_dc_tg} "
    run_cmd "${CMD1} && ${CMD2}"
  done
done


# CREATE SYSTEM BIOLOGY files --------

# for the following species...
for SPECIES in "${SPECIES_LIST[@]}"
do
  # get local variables
  OUTNAME="${SPECIES}_${DATE}"
  OUTFILE="${OUTDIR}/${OUTNAME}.categories.tsv"
  LOGFILE="${LOGDIR}/create_sb.${OUTNAME}.log"
  OUTFILE_stats="${OUTDIR}/${OUTNAME}.stats.tsv"

  OUTFILE_pid2cat="${OUTDIR}/${OUTNAME}.pid2cat.tsv"
  LOGFILE_pid2cat="${LOGDIR}/createRels.${OUTNAME}.pid2cat.log"

  # execute commands
  CMD1="python '${CODEDIR}/src/create_sb.py' -s ${SPECIES} -o '${OUTFILE}' -vv  &> '${LOGFILE}' "
  CMD2="python '${CODEDIR}/src/createRels.v0211.py' -vv  -ii '${OUTFILE}' -o '${OUTFILE_pid2cat}' -i 'Protein' -j 'cat_*' &> '${LOGFILE_pid2cat}'"
  CMD3="python '${CODEDIR}/src/stats_sb.py' -i '${OUTFILE}' -o ${OUTFILE_stats} -vv  &>> '${LOGFILE}' "
  run_cmd "${CMD1} && ${CMD2} && ${CMD3}"
done

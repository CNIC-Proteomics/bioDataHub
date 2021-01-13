#!/usr/bin/bash

# Declare variables
DATE="$(date +"%Y%m")" # with date folder
CODEDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd -P)/.."
BASEDIR="//tierra.cnic.es/SC/U_Proteomica/UNIDAD/Databases/UniProt_sync"
OUTDIR="${BASEDIR}/${DATE}" # with date folder
WSDIR="${BASEDIR}/current_release"
LOGDIR="${CODEDIR}/logs/${DATE}"

# Function that executes the input command
run_cmd () {
  echo "-- $1"  
  # eval $1
}

# prepare workspaces
mkdir "${LOGDIR}"

# for the following species...
# create the System biology database
TYPE_LIST=(sw sw-tr)
for TYPE in "${TYPE_LIST[@]}"
do
  SPECIES_LIST=(human mouse rat pig rabbit zebrafish)
  for SPECIES in "${SPECIES_LIST[@]}"
  do
    # get local variables
    LOGFILE="${LOGDIR}/create_db_sb.${SPECIES}-${TYPE}.log"
    # execute commands
    CMD="time python '${CODEDIR}/src/create_db_sb.py' -s ${SPECIES} -f ${TYPE} -o '${OUTDIR}' -vv  &> '${LOGFILE}' "
    run_cmd "${CMD}"
  done
done

# for the following species...
# create the Target/Decoy database
for INFILE in $(ls -1 "${OUTDIR}"/*.fasta)
do
    # get local variables
    filename=$(basename "${INFILE}")
    filename="${filename%.*}"
    OUTFILE_dc="${OUTDIR}/${filename}.decoy.fasta"
    OUTFILE_tg="${OUTDIR}/${filename}.target.fasta"
    OUTFILE="${OUTDIR}/${filename}.target-decoy.fasta"
    LOGFILE="${LOGDIR}/decoyPYrat.${filename}.log"
    # execute commands
    CMD="time python '${CODEDIR}/src/decoyPYrat.v2.py' --output_fasta '${OUTFILE_dc}' --decoy_prefix=DECOY '${INFILE}' &> '${LOGFILE}' && cat ${OUTFILE_tg} ${OUTFILE_dc} > ${OUTFILE}"
    run_cmd "${CMD}"
done

# for the following species...
# create the relationship file "commentline2category" (the version of protein2category)
for INFILE in $(ls -1 "${OUTDIR}"/*.categories.tsv)
do
    # get local variables
    filename=$(basename "${INFILE}")
    filename="${filename%.categories.tsv}"
    OUTFILE="${OUTDIR}/${filename}.cat.tsv"
    LOGFILE="${LOGDIR}/createRels.${filename}.log"
    # execute commands
    CMD="time python '${CODEDIR}/src/createRels.v023.py' -vv  -ii '${INFILE}' -ji '${INFILE}' -o '${OUTFILE}' -i 'Comment_Line' -j 'cat_*' &> '${LOGFILE}'"
    run_cmd "${CMD}"
done

# # Delete the last version
# mv  "${WSDIR}"  "BAKbefore_${DATE}"

# # prepare workspaces
# mkdir "${WSDIR}"

# # Copy the new version to the folder
# cp -r "${OUTDIR}/." "${WSDIR}/."


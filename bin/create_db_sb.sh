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
CTERMS=(
  "GO,cat_GO_C:cat_GO_F:cat_GO_P"
  "KEGG,cat_KEGG"
  "PANTHER,cat_PANTHER"
  "Reactome,cat_Reactome"
  "CORUM,cat_CORUM"
  "MIM,cat_OMIM"
  "DrugBank,cat_DrugBank"
  "GO_KEGG_PANTHER_Reactome,cat_GO_C:cat_GO_F:cat_GO_P:cat_KEGG:cat_PANTHER:cat_Reactome"
)


# Function that executes the input command
run_cmd () {
  echo "-- $1"
  echo ""
  eval $1
}

# prepare workspaces
mkdir "${LOGDIR}"



# CREATE FASTA files --------

# go through the repositories...
for TYPE in "${TYPE_LIST[@]}"
do
  # go through the species...
  for SPECIES in "${SPECIES_LIST[@]}"
  do
    # get local variables
    OUTDIR_spe="${OUTDIR}/${SPECIES}/sequences" # re-declare the outdir with date+version+species
    OUTNAME="${SPECIES}_${DATE}_${TYPE}"
    OUTFILE="${OUTDIR_spe}/${OUTNAME}.fasta"
    LOGFILE="${LOGDIR}/create_fasta.${OUTNAME}.log"

    OUTFILE_dc="${OUTDIR_spe}/${OUTNAME}.decoy.fasta"
    OUTFILE_tg="${OUTDIR_spe}/${OUTNAME}.target.fasta"
    OUTFILE_dc_tg="${OUTDIR_spe}/${OUTNAME}.target-decoy.fasta"
    LOGFILE_dc_tg="${LOGDIR}/decoyPYrat.${OUTNAME}.log"

    # execute commands
    CMD1="python '${CODEDIR}/src/create_fasta.py' -s ${SPECIES} -f ${TYPE} -o '${OUTFILE}' -d -vv  &> '${LOGFILE}' "
    CMD2="python '${CODEDIR}/src/decoyPYrat.v2.py' --output_fasta '${OUTFILE_dc}' --decoy_prefix=DECOY -t '${OUTFILE}.tmp' '${OUTFILE}' &> '${LOGFILE_dc_tg}' && cat ${OUTFILE_tg} ${OUTFILE_dc} > ${OUTFILE_dc_tg} "
    run_cmd "${CMD1} && ${CMD2}"
  done
done



# CREATE SYSTEM BIOLOGY file (from UniProt) and the PROTEIN2CATEGORY files --------

# go through the species...
for SPECIES in "${SPECIES_LIST[@]}"
do
  # get local variables
  OUTDIR_spe="${OUTDIR}/${SPECIES}/categories" # re-declare the outdir with date+version+species
  OUTNAME="${SPECIES}_${DATE}"
  CATFILE="${OUTDIR_spe}/${OUTNAME}.uniprot.tsv"
  STAFILE="${OUTDIR_spe}/${OUTNAME}.stats.tsv"
  LOGFILE="${LOGDIR}/create_sb.${OUTNAME}.log"

  # execute the program that creates the large file with categories from UniProt, and retrieve some statistical values.
  CMD1="python '${CODEDIR}/src/create_sb.py' -s ${SPECIES} -o '${CATFILE}' -vv  &> '${LOGFILE}' "
  CMD2="python '${CODEDIR}/src/stats_sb.py' -i '${CATFILE}' -o ${STAFILE} -vv  &>> '${LOGFILE}' "
  run_cmd "${CMD1} && ${CMD2}"

  # go through the categories...
  for CTERM in "${CTERMS[@]}"
  do IFS=","
    # get the values from the split by ','
    set ${CTERM}
    CNAME="${1}"
    CCOLS="${2}"
    # get local variables
    CNAME=$(echo ${CNAME} | tr '[:upper:]' '[:lower:]') # convert to lowercase
    OUTNAME="q2c__${SPECIES}_${DATE}.${CNAME}"
    RTFILE="${OUTDIR_spe}/${OUTNAME}.tsv"
    LOGFILE="${LOGDIR}/create_rt.${OUTNAME}.log"

    # execute the program that creates the relation table 'q2c' from the given columns of categories
    CMD3="python '${CODEDIR}/src/create_rt.py' -vv  -ii '${CATFILE}' -o '${RTFILE}' -i 'Protein' -j '${CCOLS}' &> '${LOGFILE}'"
    run_cmd "${CMD3}"
  done

done

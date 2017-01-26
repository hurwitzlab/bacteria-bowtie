#!/bin/bash

####Adjust memory as needed for size of input####

#PBS -W group_list=bhurwitz
#PBS -q standard
#PBS -l jobtype=cluster_only
#PBS -l select=1:ncpus=2:mem=3gb:pcmem=2gb
#PBS -l pvmem=6gb
#PBS -l walltime=24:00:00
#PBS -l cput=24:00:00
#PBS -M scottdaniel@email.arizona.edu
#PBS -m bea

source /usr/share/Modules/init/bash

module load R/3.2.1

set -u

COMMON="$WORKERS_DIR/common.sh"

if [ -e $COMMON ]; then
  . "$COMMON"
else
  echo Missing common \"$COMMON\"
  exit 1
fi

if [ -z $SCRIPT_DIR ]; then
  echo Missing SCRIPT_DIR
  exit 1
fi

echo Started $(date)

echo Host $(hostname)

if [[ $RUNKEGG -eq 1 ]]; then
    echo Working on Kegg...

    cd "$UPROC_ID2FEATURE/KEGG"

    TMP_FILES=$(mktemp)

    get_lines $FILES_LIST $TMP_FILES ${PBS_ARRAY_INDEX:=1} $STEP_SIZE

    NUM_FILES=$(lc $TMP_FILES)

    echo Processing \"$NUM_FILES\" input directories

    while read SAMPLE; do
        
        cd "$UPROC_ID2FEATURE/KEGG/$SAMPLE"

        Rscript --verbose $WORKERS_DIR/sumCounts.R $PWD $SAMPLE Kegg $UPROC_FEATSUM

    done < "$TMP_FILES"

fi

if [[ $RUNPFAM -eq 1 ]]; then
    echo Working on Pfam...

    cd "$UPROC_ID2FEATURE/PFAM"

    TMP_FILES=$(mktemp)

    get_lines $FILES_LIST $TMP_FILES ${PBS_ARRAY_INDEX:=1} $STEP_SIZE

    NUM_FILES=$(lc $TMP_FILES)

    echo Processing \"$NUM_FILES\" input directories

    while read SAMPLE; do
        
        cd "$UPROC_ID2FEATURE/PFAM/$SAMPLE"

        Rscript --verbose $WORKERS_DIR/sumCounts.R $PWD $SAMPLE Pfam $UPROC_FEATSUM

    done < "$TMP_FILES"

fi

echo Finished $(date)

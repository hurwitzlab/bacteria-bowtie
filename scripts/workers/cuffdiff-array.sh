#!/usr/bin/env bash

#PBS -W group_list=bhurwitz
#PBS -q qualified
#PBS -l select=1:ncpus=12:mem=36gb
###and the amount of time required to run it
#PBS -l walltime=24:00:00
#PBS -l cput=24:00:00
#PBS -M scottdaniel@email.arizona.edu
#PBS -m bea

if [ -n "$PBS_O_WORKDIR" ]; then
    cd $PBS_O_WORKDIR
fi

set -u

COMMON="$WORKER_DIR/common.sh"

if [ -e $COMMON ]; then
  . "$COMMON"
else
  echo Missing common \"$COMMON\"
  exit 1
fi

TMP_FILES=$(mktemp)

get_lines $FILES_TO_PROCESS $TMP_FILES $PBS_ARRAY_INDEX $STEP_SIZE

NUM_FILES=$(lc $TMP_FILES)

echo Found \"$NUM_FILES\" files to process

cd $ALLBACT_OUT

echo Running cuffdiff with gff $GFF in $ALLBACT_OUT

export ALLBAMS="$ALLBACT_OUT/allbams"

find $ALLBACT_OUT -iname \*.bam -print | sort > $ALLBAMS

while read GFF; do
    
    if [[ ! -d $(basename $GFF) ]]; then
        mkdir -p $(basename $GFF)
    fi

    time cuffdiff -p 12 -L S1,S2,S3,S4 \
        -o $(basename $GFF)
        --no-diff \
        --no-js-tests \
        -c 1000 \
        --max-mle-iterations 1 \
        --quiet \
        --no-update-check \
        $GFF $(cat $ALLBAMS) 

done < "$TMP_FILES"

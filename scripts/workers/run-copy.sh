#!/usr/bin/env bash

#PBS -W group_list=bhurwitz
#PBS -q standard
#PBS -l select=1:ncpus=2:mem=3gb
#PBS -l walltime=0:30:00
#PBS -l cput=0:30:00
#PBS -M scottdaniel@email.arizona.edu
#PBS -m bea

cd $PBS_O_WORKDIR

COMMON="$WORKER_DIR/common.sh"

if [ -e $COMMON ]; then
  . "$COMMON"
else
  echo Missing common \"$COMMON\"
  exit 1
fi

TMP_FILE=$(mktemp)

get_lines $STRAINS $TMP_FILE $PBS_ARRAY_INDEX $STEP_SIZE

NUM_FILES=$(lc $TMP_FILE)

echo Found \"$NUM_FILES\" files to copy

if [ -d $ALLPATRIC_ANNOT ]; then
    echo Copying files to $ALLPATRIC_ANNOT
else
    init_dir "$ALLPATRIC_ANNOT"
    echo Made $ALLPATRIC_ANNOT
fi

for FILE in $(cat $TMP_FILE); do
    set -x
    #already got the genomes
    #cp $PATRIC_GENOMES/"$FILE".fna $ALLPATRIC_ANNOT
    cp $PATRIC_ANNOT/"$FILE".PATRIC.cds.tab $ALLPATRIC_ANNOT
done

printf "\nDONE\n\n"

#!/usr/bin/env bash

#PBS -W group_list=bhurwitz
#PBS -q standard
#PBS -l jobtype=cluster_only
#PBS -l select=1:ncpus=12:mem=23gb:pcmem=2gb
#PBS -l pvmem=46gb
#PBS -l walltime=3:00:00
#PBS -l cput=36:00:00
#PBS -M scottdaniel@email.arizona.edu
#PBS -m bea

cd $PBS_O_WORKDIR

set -u

echo "Started at $(date) on host $(hostname)"

COMMON="$WORKER_DIR/common.sh"

if [ -e $COMMON ]; then
  . "$COMMON"
else
  echo Missing common \"$COMMON\"
  exit 1
fi

LEFT_TMP_FILES=$(mktemp)
RIGHT_TMP_FILES=$(mktemp)

get_lines $LEFT_FILES_LIST $LEFT_TMP_FILES $FILE_START $STEP_SIZE
get_lines $RIGHT_FILES_LIST $RIGHT_TMP_FILES $FILE_START $STEP_SIZE

NUM_FILES=$(lc $LEFT_TMP_FILES)

echo Found \"$NUM_FILES\" files to process

echo Processing $(cat $LEFT_TMP_FILES) and $(cat $RIGHT_TMP_FILES)

while read LEFT_FASTQ; do

    while read RIGHT_FASTQ; do

        test2=$(echo $RIGHT_FASTQ | sed 's/\.[1-2]\.fastq//')
        test1=$(echo $LEFT_FASTQ | sed 's/\.[1-2]\.fastq//')

        if [ "$test1" = "$test2" ]; then
            IN_LEFT=$FASTQ_DIR/$LEFT_FASTQ
            IN_RIGHT=$FASTQ_DIR/$RIGHT_FASTQ
            IN_L_UNP=$FASTQ_DIR/"$test1".nomatch1.fastq
            IN_R_UNP=$FASTQ_DIR/"$test2".nomatch2.fastq
           
            OUT_DIR=$BOWTIE2_OUT_DIR/$(dirname $LEFT_FASTQ)

            OUT=$OUT_DIR/$(basename $LEFT_FASTQ ".fa").sam

            if [[ ! -d "$OUT_DIR" ]]; then
                mkdir -p "$OUT_DIR"
            fi
            
            if [[ -e $OUT ]]; then
                echo "Sam file already exists, skipping..."
                continue
            else
                echo "Processing $LEFT_FASTQ"
            fi

            bowtie2 -p 12 \
                --very-sensitive-local \
                --no-unal \
                --no-sq \
                -k 1 \
                -x $BOWTIE2_DB \
                -1 $IN_LEFT \
                -2 $IN_RIGHT \
                -U "$IN_L_UNP","$IN_R_UNP" \
                -S $OUT
        else
            continue
        fi

        done < "$RIGHT_TMP_FILES"

done < "$LEFT_TMP_FILES"

echo "Done at $(date)"


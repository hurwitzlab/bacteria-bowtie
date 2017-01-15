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

echo "Doing sample $SAMPLE"
OUT_DIR="$MOUSE_OUT/$SAMPLE-cluster"
init_dir "$OUT_DIR"

cd "$FASTQ_DIR"

export BOWTIE2_INDEXES="$(dirname $MOUSEBT2)/"
echo "Using indexes here: $BOWTIE2_INDEXES"
echo "With basename $(basename $MOUSEBT2)"

time tophat --max-multihits 1 --read-mismatches 0 \
    -o $OUT_DIR -p 12 \
    --mate-inner-dist 25 --library-type fr-unstranded \
    --transcriptome-index=$MOUSETRANS/known --transcriptome-only \
    $(basename $MOUSEBT2) \
    $(cat $LEFT_FASTQ),$(cat $UNPAIRED) $(cat $RIGHT_FASTQ)

rm $OUT_DIR/unmapped.bam

echo "Done at $(date)"


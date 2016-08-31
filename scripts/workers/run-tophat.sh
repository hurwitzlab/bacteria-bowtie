#!/usr/bin/env bash
#PBS -m bea
#PBS -M scottdaniel@email.arizona.edu
#PBS -W group_list=bhurwitz
#PBS -q standard

### Set the number of cpus that will be used.
#PBS -l select=1:ncpus=6:mem=30Gb
#PBS -l cput=24:0:0
#PBS -l walltime=4:0:0

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
OUT_DIR="$MOUSE_OUT/$SAMPLE"
init_dir "$OUT_DIR"

cd "$FASTQ_DIR"

export BOWTIE2_INDEXES=$(dirname $MOUSEBT2)
echo "Using indexes here: $BOWTIE2_INDEXES\n \
    With basename $(basename $MOUSEBT2)"

time tophat --max-multihits 1 --no-discordant --no-mixed --read-mismatches 0 \
    -o $OUT_DIR -p 6 \
    --mate-inner-dist 25 --library-type fr-unstranded genome \
    --transcriptome-index=$MOUSETRANS/known --transcriptome-only \
    $(basename $MOUSEBT2) \
    $LEFT_FASTQ,$UNPAIRED $RIGHT_FASTQ

rm $OUT_DIR/unmapped.bam

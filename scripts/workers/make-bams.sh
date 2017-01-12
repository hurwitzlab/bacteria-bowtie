#!/usr/bin/env bash

#PBS -W group_list=bhurwitz
#PBS -q standard
#PBS -l select=1:ncpus=6:mem=6gb
#PBS -l walltime=3:00:00
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

if [ $SAMDIR == $MOUSE_OUT ]; then
    GENOME=$MOUSEFASTA
elif [ $SAMDIR == $BFRAG_OUT ]; then
    GENOME=$BFRAGFASTA
elif [ $SAMDIR == $ALLBACT_OUT ]; then
    GENOME=$ALLFASTA
elif [ $SAMDIR == $COMB_OUT ]; then
    GENOME=$COMBFASTA
fi

cd $SAMDIR

echo Converting $SAMPLE using reference $GENOME

samtools view -@ 6 -bT $GENOME $SAMPLE.sam > $SAMPLE.temp

echo Sorting $SAMPLE

samtools sort -@ 6 $SAMPLE.temp > $SAMDIR/$SAMPLE.bam

echo Removing $SAMPLE.temp

rm $SAMPLE.temp


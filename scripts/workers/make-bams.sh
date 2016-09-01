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
fi

if [ $SAMDIR == $BFRAG_OUT ]; then
    GENOME=$BFRAGFASTA
fi

echo Converting $SAMPLE using reference $GENOME

samtools view -@ 12 -bT $GENOME $SAMPLE > $SAMPLE.temp

echo Sorting $SAMPLE

samtools sort -@ 12 $SAMPLE.temp > $SAMDIR/$SAMPLE.bam

echo Removing $SAMPLE.temp

rm $SAMPLE.temp


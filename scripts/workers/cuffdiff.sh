#!/usr/bin/env bash

#PBS -W group_list=bhurwitz
#PBS -q qualified
###this is probably the memory required for the entire cuffdiff
###on a 8 million gene set
#PBS -l select=1:ncpus=12:mem=504gb:pcmem=8gb
###and the amount of time required to run it
#PBS -l walltime=180:00:00
#PBS -l cput=234:00:00
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
    GFF=$MOUSEGFF
elif [ $SAMDIR == $BFRAG_OUT ]; then
    GFF=$BFRAGGFF
elif [ $SAMDIR == $ALLBACT_OUT ]; then
    GFF=$ALLGFF
fi

cd $SAMDIR

echo Running cuffdiff with gff $GFF in $SAMDIR

export ALLBAMS="$SAMDIR/allbams"

find $SAMDIR -iname \*.bam -print | sort > $ALLBAMS

time cuffdiff -p 12 -L S1,S2,S3,S4 $GFF $(cat $ALLBAMS) 

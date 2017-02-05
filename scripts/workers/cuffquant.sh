#!/usr/bin/env bash

#PBS -W group_list=bhurwitz
#PBS -q qualified
#PBS -l select=1:ncpus=12:mem=72gb
###and the amount of time required to run it
#PBS -l walltime=24:00:00
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

set -x

if [ $SAMDIR == $MOUSE_OUT ]; then
    GFF=$MOUSEGFF
elif [ $SAMDIR == $BFRAG_OUT ]; then
    GFF=$BFRAGGFF
elif [ $SAMDIR == $ALLBACT_OUT ]; then
    export GFF=$ALLGFF
    export rRNAGFF=$ALLrRNAGFF
elif [ $SAMDIR == $COMB_OUT ]; then
    export GFF=$COMBGFF
    export rRNAGFF=$COMBrRNAGFF
elif [ $SAMDIR == $UNK_OUT ]; then
    export GFF=$UNKGFF
    export rRNAGFF=$UNKrRNAGFF
fi

cd $SAMDIR

echo Running cuffdiff with gff $GFF in $SAMDIR
echo With sample "$SAMPLE".bam

#export ALLBAMS="$SAMDIR/allbams"

#find $SAMDIR -iname \*.bam -print | sort > $ALLBAMS

if [[ ! -d "$SAMPLE"-cuffquant-out ]]; then
    mkdir -p "$SAMPLE"-cuffquant-out
else
    rm "$SAMPLE"-cuffquant-out/*
fi

time cuffquant -p 12 \
    -o "$SAMPLE"-cuffquant-out \
    -M $rRNAGFF \
    --quiet \
    --no-length-correction \
    $GFF "$SAMPLE".bam

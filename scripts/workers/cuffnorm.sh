#!/usr/bin/env bash

#PBS -W group_list=bhurwitz
#PBS -q standard
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

export ALLCXBS="$SAMDIR/allcxbs"

find $SAMDIR -iname \*.cxb -print | sort > $ALLCXBS

if [[ ! -d cuffnorm-out ]]; then
    mkdir -p cuffnorm-out
else
    rm -r cuffnorm-out/*
fi

time cuffnorm -p 12 --labels S1,S2,S3,S4 \
    -o combined-cuffnorm-out \
    --quiet \
    $GFF $(cat $ALLCXBS)

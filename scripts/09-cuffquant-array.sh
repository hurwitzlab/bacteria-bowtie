#!/usr/bin/env bash

#
# This script is intended to do cuffquant for very large genomes/annotation files by splitting up the annotation file into chunks
#
unset module
set -u
source ./config.sh
export CWD="$PWD"
export STEP_SIZE=1
#SPLIT_SIZE is number of lines to split each gff
export SPLIT_SIZE=20
echo Setting up log files...

PROG=`basename $0 ".sh"`
#Just going to put stdout and stderr together into stdout
STDOUT_DIR="$CWD/out/$PROG"

init_dir "$STDOUT_DIR"

export GFFLIST="$PRJ_DIR/tmp/gfflist"

find $(dirname $ALLGFF) -iname \*splitgff\* > $GFFLIST

export NUM_FILES=$(lc $GFFLIST)

if [[ $NUM_FILES -eq 0 ]]; then
    cd $(dirname $ALLGFF)

    #PROTIP: l/$STEP_SIZE tells split to NOT split within lines!
    split -e -d -n l/$SPLIT_SIZE $(basename $ALLGFF) splitgff
        
    find $(dirname $ALLGFF) -iname \*splitgff\* > $GFFLIST

    export NUM_FILES=$(lc $GFFLIST)

fi

echo Doing these many gffs: $NUM_FILES

qsub -J 1-$NUM_FILES:$STEP_SIZE -V -j oe -o "$STDOUT_DIR" $WORKER_DIR/cuffquant-array.sh

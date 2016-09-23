#!/usr/bin/env bash

#
# This script is intended to make bams from sams and then sort
#
unset module
set -u
source ./config.sh
export CWD="$PWD"
export STEP_SIZE=1
echo Setting up log files...

PROG=`basename $0 ".sh"`
#Just going to put stdout and stderr together into stdout
STDOUT_DIR="$CWD/out/$PROG"

init_dir "$STDOUT_DIR"

export GFFLIST="$PRJ_DIR/tmp/gfflist"

find $(dirname $ALLGFF) -iname \*splitgff\* > $GFFLIST

export NUM_FILES=$(lc $GFFLIST)

echo Doing these many gffs: $NUM_FILES

qsub -J 1-$NUM_FILES:$STEP_SIZE -V -j oe -o "$STDOUT_DIR" $WORKER_DIR/cuffdiff-array.sh

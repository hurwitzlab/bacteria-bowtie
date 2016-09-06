#!/usr/bin/env bash

#
# This script is intended to make bams from sams and then sort
#
#echo "Number of arguments is $#"
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

export SAMDIR=$MOUSE_OUT

qsub -V -j oe -o "$STDOUT_DIR" $WORKER_DIR/cuffdiff.sh

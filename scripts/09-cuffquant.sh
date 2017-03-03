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

if [[ $# = 0 ]]; then
    echo "Need to know what sam output directory to work on"
    echo "e.g. ./07-cuffquant.sh \$MOUSE_OUT"
    echo "ALSO: did you source config.sh?"
    exit 1
fi

echo Setting up log files...
PROG=`basename $0 ".sh"`
#Just going to put stdout and stderr together into stdout
STDOUT_DIR="$CWD/out/$PROG"

init_dir "$STDOUT_DIR/$(basename $1)"

echo Submitting job...

for i in $SAMPLE_NAMES; do
    export SAMDIR=$1
    export SAMPLE=$i
    echo $i
    qsub -V -j oe -o "$STDOUT_DIR/$(basename $1)" $WORKER_DIR/cuffquant.sh
done


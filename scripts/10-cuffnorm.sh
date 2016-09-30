#!/usr/bin/env bash

#
# This script is intended to run cuffnorm
#
unset module
set -u
source ./config.sh
export CWD="$PWD"

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

export SAMDIR=$1

qsub -V -j oe -o "$STDOUT_DIR/$(basename $1)" $WORKER_DIR/cuffnorm.sh

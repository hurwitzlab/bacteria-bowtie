#!/usr/bin/env bash
#
# Intended to make a transcriptome index first
#

unset module
source ./config.sh

CWD=$(pwd)
PROG=`basename $0 ".sh"`
STDOUT_DIR="$CWD/out/$PROG"

init_dir "$STDOUT_DIR" "$MOUSETRANS"

JOB=$(qsub -V -N tophattransm -j oe -o "$STDOUT_DIR" $WORKER_DIR/run-transcrpt.sh)

if [ $? -eq 0 ]; then
    echo Submitted job \"$JOB\" for you. Weeeeeeeeeeeee!
else
    echo -e "\nError submitting job\n$JOB\n"
fi


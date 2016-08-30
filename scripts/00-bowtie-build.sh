#!/usr/bin/env bash
# 
# Intended to build a bowtie2 index for you
#

#all this does is launch bowite2-build for you
unset module
set -u
export CWD="$PWD"
CONFIG="./config.sh"

if [ -e $CONFIG ]; then
    . "$CONFIG"
else
    echo MIssing config \"$CONFIG\"
    exit 12385
fi

PROG=`basename $0 ".sh"`
#Just going to put stdout and stderr together into stdout
STDOUT_DIR="$CWD/out/$PROG"

init_dir "$STDOUT_DIR"


COMMON="$WORKER_DIR/common.sh"

if [ -e $COMMON ]; then
  . "$COMMON"
else
  echo Missing common \"$COMMON\"
  exit 1
fi

JOB=$(qsub -V -N bowtie2-build -j oe -o "$STDOUT_DIR" $WORKER_DIR/run-bowtie2-build.sh)

if [ $? -eq 0 ]; then
    echo Submitted job \"$JOB\" for you. Weeeeeeeeeeeee!
else
    echo -e "\nError submitting job\n$JOB\n"
fi


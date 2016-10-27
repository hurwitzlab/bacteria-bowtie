#!/usr/bin/env bash

#
#This script is intended to grab those muy importante contigs
#
unset module
set -u
source ./config.sh
export CWD="$PWD"
export STEP_SIZE=100
#lines of just-strain-list.txt= 1944
#divide by 100 and you get 20 subjobs

echo Setting up log files...
PROG=`basename $0 ".sh"`
#Just going to put stdout and stderr together into stdout
export STDOUT_DIR="$CWD/out/$PROG"

init_dir "$STDOUT_DIR"

if [ -e $STRAINS ]; then
    echo "mapping file $STRAINS exists, assuming it's fine"
else
    echo "no $STRAINS or bad file, exiting..."
    exit 12345
fi

export NUM_COPIES=$(lc $STRAINS)

echo \"Processing $NUM_COPIES contigs\"

JOB=$(qsub -J 1-$NUM_COPIES:$STEP_SIZE -V -N fetch_dog -j oe -o "$STDOUT_DIR" $WORKER_DIR/run-copy.sh)

if [ $? -eq 0 ]; then
  echo Submitted job \"$JOB\" for you. La la la. La. la.
else
  echo -e "\nError submitting job\n$JOB\n"
fi





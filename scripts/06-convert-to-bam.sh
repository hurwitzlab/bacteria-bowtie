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
    echo "e.g. ./05-convert-to-bam.sh \$MOUSE_OUT"
    echo "ALSO: did you source config.sh?"
    exit 1
fi

echo Setting up log files...
PROG=`basename $0 ".sh"`
#Just going to put stdout and stderr together into stdout
STDOUT_DIR="$CWD/out/$PROG"

init_dir "$STDOUT_DIR/$(basename $1)"
#
#if [ $1 == $MOUSE_OUT ]; then
#    echo "You have chosen Mouse"
#    echo "Directory is $MOUSE_OUT"
#    echo "Files there are $(ls $MOUSE_OUT)"
##    exit 1
#fi
#
#if [ $1 == $BFRAG_OUT ]; then
#    echo "You have chosen Bfrag"
#    echo "Directory is $BFRAG_OUT"
#    echo "Files there are $(ls $BFRAG_OUT)"
#fi
#
#export FILE_LIST="$TEMP_DIR/sam_files"
#
#find $1 -iname \*RNA\*sam > $FILE_LIST
#
#export NUM_FILES=$(lc $FILE_LIST)
#
#echo \"Found $NUM_FILES to make bams from\"
#echo \"Splitting them up in batches of "$STEP_SIZE"\"
echo Submitting job...
 
for i in $SAMPLE_NAMES; do
    export SAMDIR=$1
    export SAMPLE=$i
    echo $1
    echo $i
    qsub -V -j oe -o "$STDOUT_DIR/$(basename $1)" $WORKER_DIR/make-bams.sh
done

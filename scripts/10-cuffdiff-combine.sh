#!/usr/bin/env bash

#
# This script is intended to combine cuffdiff outputs
#
unset module
set -u
source ./config.sh
export CWD="$PWD"

echo Setting up log files...

PROG=`basename $0 ".sh"`
#Just going to put stdout and stderr together into stdout
STDOUT_DIR="$CWD/out/$PROG"

init_dir "$STDOUT_DIR"

cd $ALLBACT_OUT

export GENECOUNTS="isoforms.fpkm_tracking"

if [[ $(wc -l $(find ./ -iname \*$GENECOUNTS)) -ge 1 ]]; then
    find ./ -iname \*$GENECOUNTS -print0 | xargs -0 -I file cat file \
        > total.fpkm_tracking
else
    echo "YOu got the name wrong"
    echo "Or couldn't find anything"
    exit 1
fi




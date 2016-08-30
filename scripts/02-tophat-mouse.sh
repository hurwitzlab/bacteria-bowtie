#
# to map the mouse
#

source ./config.sh

CWD=$(pwd)
PROG=`basename $0 ".sh"`
STDOUT_DIR="$CWD/out/$PROG"

init_dir "$STDOUT_DIR" "$MOUSE_OUT"

cd "$FASTQ_DIR"

for i in $SAMPLE_NAMES; do

    export LEFT_FASTQ="$TEMP_DIR/$i-left_fastq"
    export RIGHT_FASTQ="$TEMP_DIR/$i-right-fastq"
    export UNPAIRED="$TEMP_DIR/$i-unpaired-fastq"

    find . -type f -regextype 'sed' -iregex "\.\/$i.*\.1.fastq" | sort > $LEFT_FASTQ
    find . -type f -regextype 'sed' -iregex "\.\/$i.*\.2.fastq" | sort > $RIGHT_FASTQ
    find . -type f -regextype 'sed' -iregex "\.\/$i.*nomatch.*" > $UNPAIRED

    echo "Mapping $i FASTQs to $MOUSEBT2"

done

JOB=$(qsub -V -N tophatm -j oe -o "$STDOUT_DIR" $WORKER_DIR/tophat.sh)

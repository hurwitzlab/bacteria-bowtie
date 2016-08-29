#
# to map the mouse
#

source ./config.sh

PROG=`basename $0 ".sh"`
STDOUT_DIR="$CWD/out/$PROG"

init_dir "$STDOUT_DIR" "$MOUSE_OUT"

cd "$FASTQ_DIR"

export LEFT_FASTQ="$PRJ_DIR/left_fastq"
export RIGHT_FASTQ="$PRJ_DIR/right_fastq"
export UNPAIRED="$PRJ_DIR/unpaired_fastq"

find . -type f -regextype 'sed' -iregex '.+1.fastq.+' > $LEFT_FASTQ
find . -type f -regextype 'sed' -iregex '.+2.fastq.+' > $RIGHT_FASTQ
find . -type f -regextype 'sed' -iregex '.+nomatch.+' > $UNPAIRED

echo "Mapping FASTQs to $MOUSEBT2"

JOB=$(qsub -V -N tophatm -j oe -o "$STDOUT_DIR" $WORKER_DIR/tophat.sh)

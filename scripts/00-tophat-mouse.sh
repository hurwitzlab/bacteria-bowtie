#
# This script is intended to get back the original fastq's from your screened fasta's
#

source ./config.sh
export STEP_SIZE=20

PROG=`basename $0 ".sh"`
STDOUT_DIR="$CWD/out/$PROG"

init_dir "$STDOUT_DIR" "$FILTERED_FQ"

cd "$FASTA_DIR"

export FILES_LIST="$FASTA_DIR/files-list"

find . -type f -name \*.fa | sed "s/^\.\///" > $FILES_LIST

NUM_FILES=$(lc $FILES_LIST)

echo Found \"$NUM_FILES\" files in \"$FASTA_DIR\"

JOB=$(qsub -J 1-$NUM_FILES:$STEP_SIZE -V -N fetch-fq -j oe -o "$STDOUT_DIR" $WORKER_DIR/fetch-fq.sh)

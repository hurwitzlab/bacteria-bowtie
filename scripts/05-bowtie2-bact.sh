#
# to map the bacteria
#
unset module
source ./config.sh

CWD=$(pwd)
PROG=`basename $0 ".sh"`
STDOUT_DIR="$CWD/out/$PROG"

init_dir "$STDOUT_DIR" "$COMBINED_OUT"

cd "$FASTQ_DIR"

for i in $SAMPLE_NAMES; do
    
    export SAMPLE=$i
    export LEFT_FASTQ="$TEMP_DIR/$i-left-fastq"
    export RIGHT_FASTQ="$TEMP_DIR/$i-right-fastq"
    export UNPAIRED="$TEMP_DIR/$i-unpaired-fastq"
    export BT2="$COMBBT2"
    export OUT="$COMBINED_OUT"

#DON'T HAVE TO DO THIS EVERYTIME IF FILES ALREADY THERE
    find . -type f -regextype 'sed' -iregex "\.\/$i.*\.1.fastq" \
         > $LEFT_FASTQ
    
    find . -type f -regextype 'sed' -iregex "\.\/$i.*\.2.fastq"  \
         > $RIGHT_FASTQ
    
    find . -type f -regextype 'sed' -iregex "\.\/$i.*nomatch.*" \
         > $UNPAIRED 
    
    echo "Mapping $i FASTQs to $BT2"

    JOB=$(qsub -V -N bowtie3 -j oe -o "$STDOUT_DIR" $WORKER_DIR/run-bowtie2.sh)

    if [ $? -eq 0 ]; then
        echo Submitted job \"$JOB\" for you. Weeeeeeeeeeeee!
    else
        echo -e "\nError submitting job\n$JOB\n"
    fi

done


#!/bin/bash
#PBS -m bea
#PBS -M scottdaniel@email.arizona.edu
#PBS -W group_list=bhurwitz
#PBS -q standard
#PBS -l jobtype=smp_only

### Set the number of cpus that will be used.
#PBS -l select=1:ncpus=12:mem=30Gb
#PBS -l cput=32:0:0
#PBS -l walltime=2:0:0

source /usr/share/Modules/init/bash

module load bowtie2
module load tophat
module load samtools
module load cufflinks

SAMPLE=$i
OUT_DIR="$MOUSE_OUT/$SAMPLE"
init_dir "$OUT_DIR"
time tophat --max-multihits 1 --no-discordant --no-mixed --read-mismatches 0 \
    -o $OUT_DIR -p 12 \
    --mate-inner-dist 25 --library-type fr-unstranded genome \
    --transcriptome-index=$MOUSETRANS/known --transcriptome-only \
    $MOUSEBT2
    $LEFT_FASTQ,$UNPAIRED $RIGHT_FASTQ

rm $OUT_DIR/unmapped.bam

#!/bin/bash

# --------------------------------------------------
#
# config.sh
#
# Edit this file to match your directory structure
#
# --------------------------------------------------

#
# Some constants
#
export QSTAT="/usr/local/bin/qstat_local"
export GUNZIP="/bin/gunzip"

#
# The main checkout
#
export PRJ_DIR="/rsgrps/bhurwitz/scottdaniel/tophat-bacteria"

# 
# Sample names
#
export SAMPLE_NAMES="RNA_1 RNA_2 RNA_3 RNA_4"

#
# Where mouse bowtie indices are and wholegenomefasta
#
export REF_DIR="/rsgrps/bhurwitz/hurwitzlab/data/reference"
export MOUSEFASTA="$REF_DIR/tophat_mouse/Mus_musculus/NCBI/build37.2/Sequence/WholeGenomeFasta/genome.fa"
export MOUSEBT2="$REF_DIR/tophat_mouse/Mus_musculus/NCBI/build37.2/Sequence/Bowtie2Index/genome"
export MOUSEGFF="$REF_DIR/tophat_mouse/Mus_musculus/NCBI/build37.2/Annotation/Genes/genes.gtf"
export MOUSETRANS="$REF_DIR/tophat_mouse/Mus_musculus/NCBI/build37.2/transcriptome"

#
# Where bacteria genomes are
#
export BACT_DIR="$PRJ_DIR/genomes"
#bfrag
export BFRAGBT2="$BACT_DIR/bfrag/genome"
export BFRAGGFF="$BACT_DIR/bfrag/genome.gff"

#
# Where we can find the worker scripts
#
export SCRIPT_DIR="$PRJ_DIR/scripts"
export WORKER_DIR="$SCRIPT_DIR/workers"
#
# Where to put all our generated data
#
export DATA_DIR="$PRJ_DIR/data"
export MOUSE_OUT="$DATA_DIR/topmouse"
export BFRAG_OUT="$DATA_DIR/bfrag_out"

# Where our reads are (Qc'd and sorted into paired and unpaired)
export FASTQ_DIR="/rsgrps/bhurwitz/scottdaniel/mouseRNA/data/fastq-sorted"

# place to put temp stuff like lists of files
export TEMP_DIR="$PRJ_DIR/tmp"


#
# Some custom functions for our scripts
#
# --------------------------------------------------
function init_dir {
    for dir in $*; do
        if [ -d "$dir" ]; then
            rm -rf $dir/*
        else
            mkdir -p "$dir"
        fi
    done
}

# --------------------------------------------------
function lc() {
    wc -l $1 | cut -d ' ' -f 1
}

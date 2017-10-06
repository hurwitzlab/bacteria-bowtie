#!/bin/bash

# --------------------------------------------------
#
# config.sh
#
# Edit this file to match your directory structure
#
# --------------------------------------------------
#
# The main checkout
#
export PRJ_DIR="/rsgrps/bhurwitz/scottdaniel/bacteria-tcga"

# 
# Sample names
#
export SAMPLE_NAMES=""


export REF_DIR="/rsgrps/bhurwitz/hurwitzlab/data/reference"

# Where bacteria genomes are
# and annotation
#
export PAT_GENOMES="$REF_DIR"/patric_bacteria
#complete patric genomes
export PAT_GENOMES_COMP="$PAT_GENOMES"/complete
#whole genome shotgun patric genomes
export PAT_GENOMES_WGS="$PAT_GENOMES"/wgs

#patric bowtie2 indices
export PATRIC_BT2="$REF_DIR"/patric_bt2
#the complete.fa
#just using complete right now
export FA="$PATRIC_BT2"/complete.fa



#
# Where we can find the worker scripts
#
export SCRIPT_DIR="$PRJ_DIR/scripts"
export WORKER_DIR="$SCRIPT_DIR/workers"
#
# Where to put all our generated data
#
export DATA_DIR="$PRJ_DIR/data"
# place to put temp stuff like lists of files
export TEMP_DIR="$PRJ_DIR/tmp"

#where patric annotation is
export PATRIC_ANNOT=""$REF_DIR"/patric_annot"
export PAT_GFF="$PATRIC_ANNOT"/gff
export PAT_TAB="$PATRIC_ANNOT"/cdsTab

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

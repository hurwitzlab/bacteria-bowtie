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
export PRJ_DIR="/rsgrps/bhurwitz/scottdaniel/bacteria-bowtie"

# 
# Sample names
#
export SAMPLE_NAMES="RNA_1 RNA_2 RNA_3 RNA_4"

#
# Where mouse bowtie indices are and wholegenomefasta
#
export REF_DIR="/rsgrps/bhurwitz/hurwitzlab/data/reference"

#using ensembl for tophat because i just want to try it again
export MOUSEFASTA="$REF_DIR/tophat_mouse/ensembl/bowtie2/genome.fa"
export MOUSEBT2="$REF_DIR/tophat_mouse/ensembl/bowtie2/genome"
export MOUSEGFF="$REF_DIR/tophat_mouse/ensembl/bowtie2/genome.gtf"
export MOUSETRANS="$REF_DIR/tophat_mouse/ensembl/bowtie2/transcriptome"

#export MOUSEFASTA="$REF_DIR/tophat_mouse/Mus_musculus/NCBI/build37.2/Sequence/WholeGenomeFasta/genome.fa"
#export MOUSEBT2="$REF_DIR/tophat_mouse/Mus_musculus/NCBI/build37.2/Sequence/Bowtie2Index/genome"
#export MOUSEGFF="$REF_DIR/tophat_mouse/Mus_musculus/NCBI/build37.2/Annotation/Genes/genes.gtf"
#export MOUSETRANS="$REF_DIR/tophat_mouse/Mus_musculus/NCBI/build37.2/transcriptome"
#
#
# Where bacteria genomes are
# and annotation
#
export BACT_DIR="$PRJ_DIR/genomes"
#bfrag
export BFRAGBT2="$BACT_DIR/bfrag/genome"
export BFRAGGFF="$BACT_DIR/bfrag/genome.gff"
export BFRAGFASTA="$BACT_DIR/bfrag/genome.fa"
#all 1944 best strains 
export ALLBT2="$BACT_DIR/all/bowtie2_index/genome"
#export ALLGFF="$BACT_DIR/all/bowtie2_index/all.gff"
export ALLGFF="$BACT_DIR/all/bowtie2_index/all-refseq-CDS.gff"
export ALLrRNAGFF="$BACT_DIR/all/bowtie2_index/all-refseq-rRNA.gff"
export ALLFASTA="$BACT_DIR/all/bowtie2_index/all.fa"
#all those patric CDS tabs for the 1944
export ALLPATRIC_ANNOT="$BACT_DIR/all/patric_annot"

#the combined: the 1944 "best" genomes and the "low alignment score / unknown" that was assembled and then annotated
export COMB_DIR="$BACT_DIR/combined"
export COMBFASTA="$COMB_DIR/combined.fa"
export COMBBT2="$COMB_DIR/genome"
export COMBGFF="$COMB_DIR/combined.gff"
export COMBrRNAGFF="$COMB_DIR/rRNA-combined.gff"

#
# Where we can find the worker scripts
#
export SCRIPT_DIR="$PRJ_DIR/scripts"
export WORKER_DIR="$SCRIPT_DIR/workers"
#
# Where to put all our generated data
#
export DATA_DIR="$PRJ_DIR/data"
#export MOUSE_OUT="$DATA_DIR/mouse_out"
export MOUSE_OUT="$DATA_DIR/topmouse"
export BFRAG_OUT="$DATA_DIR/bfrag_out"
#export ALLBACT_OUT="$DATA_DIR/theusualsuspects"
export ALLBACT_OUT="$DATA_DIR/allbact_out"
export COMB_OUT="$DATA_DIR/combined_out"

# Where our reads are (Qc'd and sorted into paired and unpaired)
export FASTQ_DIR="/rsgrps/bhurwitz/scottdaniel/RNA_qc/data/fastq_out"

# place to put temp stuff like lists of files
export TEMP_DIR="$PRJ_DIR/tmp"

#where patric annotation is
export PATRIC_ANNOT="/rsgrps/bhurwitz/hurwitzlab/data/reference/patric_annot/patric_cds"
#where patric gff's are
export PATRIC_GFFS="/rsgrps/bhurwitz/hurwitzlab/data/reference/patric_annot/gff"
export REFSEQ_GFFS=$PATRIC_GFFS
#source file with patric accn's (first field) [space] and patric strain/genome number (second field)
export SOURCE_MAP="$DATA_DIR/just-strain-and-accn-right-list.txt"
#just the unique strains/genomes
export STRAINS="$DATA_DIR/just-strain-list.txt"



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

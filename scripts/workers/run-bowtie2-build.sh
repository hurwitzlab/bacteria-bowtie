#!/usr/bin/env bash

#PBS -W group_list=bhurwitz
#PBS -q standard
#PBS -l select=1:ncpus=12:mem=68gb
#PBS -l walltime=72:00:00
#PBS -l cput=288:00:00
#PBS -M scottdaniel@email.arizona.edu
#PBS -m bea

#make sure this matches ncpus in the above header!
export THREADS="--threads 12"
#
# runs bowtie2-build 
#

# --------------------------------------------------
# singularity is needed to run singularity images
module load singularity
# --------------------------------------------------


unset module
set -u

COMMON="$WORKER_DIR/common.sh"

if [ -e $COMMON ]; then
  . "$COMMON"
else
  echo Missing common \"$COMMON\"
  exit 1
fi

echo "Started at $(date) on host $(hostname)"

echo "Bowtie2 indexing..."

export WD="$(dirname $BT2)"
export FANAME="$(basename $FASTA)"
export BT2BASE="$(basename $BT2)"

export bt2build="singularity exec \
    -B $WD:$SING_WD \
    $SING_IMG/bowcuff.img bowtie2-build" 

cd $WD

$bt2build $THREADS --large-index $SING_WD/$FASTA $SIGN_WD/$BT2BASE
#
#if [[ ! -e "$BT2".fa ]]; then
#    ln $FASTA "$BT2".fa
#fi
#
echo "Done $(date)"

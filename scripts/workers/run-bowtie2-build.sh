#!/usr/bin/env bash

#PBS -W group_list=bhurwitz
#PBS -q standard
#PBS -l select=1:ncpus=12:mem=23gb
#PBS -l walltime=24:00:00
#PBS -l cput=24:00:00
#PBS -M scottdaniel@email.arizona.edu
#PBS -m bea

cd $PBS_O_WORKDIR

set -u

echo "Started at $(date) on host $(hostname)"

echo "Bowtie2 indexing..."

cd $(dirname $BT2)

bowtie2-build --large-index $FASTA $BT2
#
#if [[ ! -e "$BT2".fa ]]; then
#    ln $FASTA "$BT2".fa
#fi
#
echo "Done $(date)"

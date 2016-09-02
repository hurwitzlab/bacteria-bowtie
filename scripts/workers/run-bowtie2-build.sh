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

cd $(dirname $MOUSEBT2)

bowtie2-build $MOUSEFASTA $MOUSEBT2

if [[ ! -e "$MOUSEBT2".fa ]]; then
    ln $MOUSEFASTA "$MOUSEBT2".fa
fi

echo "Done $(date)"

#!/usr/bin/env bash

#PBS -W group_list=bhurwitz
#PBS -q standard
#PBS -l jobtype=cluster_only
#PBS -l select=1:ncpus=12:mem=23gb:pcmem=2gb
#PBS -l pvmem=46gb
#PBS -l walltime=24:00:00
#PBS -l cput=24:00:00
#PBS -M scottdaniel@email.arizona.edu
#PBS -m bea

cd $PBS_O_WORKDIR

set -u

echo "Started at $(date) on host $(hostname)"

echo "Bowtie2 indexing..."

cd $(dirname $MOUSEBT2)

bowtie2-build --large-index $MOUSEFASTA $MOUSEBT2

ln $MOUSEFASTA "$MOUSEBT2".fa

echo "Done $(date)"

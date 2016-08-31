#!/bin/bash
#PBS -m bea
#PBS -M scottdaniel@email.arizona.edu
#PBS -W group_list=bhurwitz
#PBS -q standard

### Set the number of cpus that will be used.
#PBS -l select=1:ncpus=16:mem=64gb
#PBS -l walltime=12:00:00
#PBS -l cput=12:00:00

echo "Time started $(date)"

time tophat  -p 16 \
    -G $MOUSEGFF \
    --transcriptome-index=$MOUSETRANS/known \
    $MOUSEBT2

echo "Time ended $(date)"

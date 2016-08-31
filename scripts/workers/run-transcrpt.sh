#!/bin/bash
#PBS -m bea
#PBS -M scottdaniel@email.arizona.edu
#PBS -W group_list=bhurwitz
#PBS -q standard

### Set the number of cpus that will be used.
#PBS -l select=1:ncpus=6:mem=1gb
#PBS -l place=free:shared
#PBS -l walltime=00:30:00

echo "Time started $(date)"

time tophat  -p 6 \
    -G $MOUSEGFF \
    --transcriptome-index=$MOUSETRANS/known \
    $MOUSEBT2

echo "Time ended $(date)"

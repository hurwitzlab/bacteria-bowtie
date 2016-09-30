#!/usr/bin/env bash

#PBS -W group_list=bhurwitz
#PBS -q qualified
#PBS -l select=1:ncpus=6:mem=36gb
#PBS -l walltime=24:00:00
#PBS -l cput=24:00:00
#PBS -M scottdaniel@email.arizona.edu
#PBS -m bea

set -u

cd /rsgrps/bhurwitz/scottdaniel/tophat-bacteria/data/allbact_out

export GFF="/rsgrps/bhurwitz/scottdaniel/tophat-bacteria/genomes/all/bowtie2_index/LPS-refseq-CDS.gff"

echo Doing RNA_1
htseq-count -a 0 -q -f bam -t CDS -i ID -s no RNA_1.bam $GFF > RNA_1.rawcounts
echo Doing RNA_2
htseq-count -a 0 -q -f bam -t CDS -i ID -s no RNA_2.bam $GFF > RNA_2.rawcounts
echo Doing RNA_3
htseq-count -a 0 -q -f bam -t CDS -i ID -s no RNA_3.bam $GFF > RNA_3.rawcounts
echo Doing RNA_4
htseq-count -a 0 -q -f bam -t CDS -i ID -s no RNA_4.bam $GFF > RNA_4.rawcounts

#!/bin/bash

#SBATCH -J transfer           # Job name
#SBATCH -o ./out/transfer.%j.out    # Specify stdout output file (%j expands to jobId)
#SBATCH -p normal           # Queue name
#SBATCH -N 1                    # Total number of nodes requested (16 cores/node)
#SBATCH -n 16                     # Total number of tasks
#SBATCH -t 24:00:00              # Run time (hh:mm:ss) - 1.5 hours
#SBATCH --mail-user=scottdaniel@email.arizona.edu
#SBATCH --mail-type=all
#SBATCH -A iPlant-Collabs         # Specify allocation to charge against

#automagic offloading for the xeon phi co-processor
#in case anything uses Intel's Math Kernel Library
export MKL_MIC_ENABLE=1
export OMP_NUM_THREADS=16
export MIC_OMP_NUM_THREADS=240
export OFFLOAD_REPORT=2

CONFIG="./config.sh"

if [[ ! -e $CONFIG ]]; then
    echo $CONFIG missing
    exit 1
else
    source $CONFIG
fi

NEW_PRJ_DIR="/work/03859/sdaniel/tophat-bacteria"

if [[ ! -d $NEW_PRJ_DIR/data/allbact_out ]]; then
    mkdir -p $NEW_PRJ_DIR/data/allbact_out
fi

if [[ ! -d $NEW_PRJ_DIR/genomes/all/bowtie2_index ]]; then
    mkdir -p $NEW_PRJ_DIR/genomes/all/bowtie2_index
fi


while [ 1 ]; do
    rsync --timeout=60 --log-file=./out/rsync-log.txt -aivvzhP --rsh=ssh --exclude='.sam' scottdaniel@login.hpc.arizona.edu:$ALLBACT_OUT $NEW_PRJ_DIR/data/allbact_out
    if [ "$?" = "0" ]; then
        echo "Rsync success!"
        exit 0
    else
        echo "Rsync failure. Retrying in a minute"
        sleep 60
    fi
done

while [ 1 ]; do
    rsync --timeout=60 --log-file=./out/rsync-log2.txt -aivvzhP --rsh=ssh --exclude='*' --include='.gff' scottdaniel@login.hpc.arizona.edu:$(dirname $ALLGFF) $NEW_PRJ_DIR/genomes/all/bowtie2_index
    if [ "$?" = "0" ]; then
        echo "Rsync success!"
        exit 0
    else
        echo "Rsync failure. Retrying in a minute"
        sleep 60
    fi
done


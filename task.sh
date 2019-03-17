#!/bin/bash

#SBATCH --job-name overlap
#SBATCH --account=dpicls
#SBATCH --output=/gscratch/ljin1/data/twitter/log/overlap.%j.out
#SBATCH --error=/gscratch/ljin1/data/twitter/log/overlap.%j.err
#SBATCH --mail-type=END,FAIL
### #SBATCH --mail-type=NONE
#SBATCH --mail-user=ljin1@uwyo.edu
#SBATCH --time=6-23:59:59

# Start the job
if [[ $OSTYPE == darwin* ]] # MacBook Pro @ libaooutrage (macOS)
then
    DATA_DIR=/Users/libao/Document/data/twitter/result
elif [[ $OSTYPE == linux-gnu ]] # Teton @ UWyo ARCC (Linux)
then
    DATA_DIR=/gscratch/ljin1/data/twitter/result
fi

# Set up output directory
if [ ! -d $DATA_DIR ]; then
    mkdir -p $DATA_DIR
fi

echo $NCPU
mpirun -n $NCPU python3 ./build_overlap_matrix.py

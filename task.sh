#!/bin/bash

#SBATCH --job-name overlap
#SBATCH --account=dpicls
#SBATCH --output=/gscratch/ljin1/data/twitter/log/overlap.%j.out
#SBATCH --error=/gscratch/ljin1/data/twitter/log/overlap.%j.err
#SBATCH --mail-type=END,FAIL
## SBATCH --mail-type=NONE
#SBATCH --mail-user=ljin1@uwyo.edu
#SBATCH --time=6-23:59:59

# Start the job
if [[ $OSTYPE == darwin* ]] # MacBook Pro @ libaooutrage (macOS)
then
    DATA_DIR=/Users/libao/Documents/data/twitter/csv
elif [[ $OSTYPE == linux-gnu ]] # Teton @ UWyo ARCC (Linux)
then
    DATA_DIR=/gscratch/ljin1/data/twitter/csv
fi

# Set up output directory
if [ ! -d $DATA_DIR ]; then
    mkdir -p $DATA_DIR
fi

echo $NCPU
module purge -q
module load arcc/0.1 slurm/18.08 swset/2018.05
if [[ $MPI_VENDOR == "intel" ]]; then
    echo "Intel MPI"
    module load intel/18.0.1 intel-mpi/2018.2.199 intel-mkl/2018.2.199
    srun python3 ./build_overlap_matrix.py -l $DATA_DIR
    #srun python3 ./build_overlap_matrix.py -l $DATA_DIR -r 1000
else
    echo "Open MPI"
    module load gcc/7.3.0 openmpi/3.1.0 intel-mkl/2018.2.199
    mpirun -np $NCPU python3 ./build_overlap_matrix.py -l $DATA_DIR
    #mpirun -np $NCPU python3 ./build_overlap_matrix.py -l $DATA_DIR -r 1000
fi

which python3

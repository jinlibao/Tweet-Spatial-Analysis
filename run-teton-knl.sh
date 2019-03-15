#!/bin/bash

#SBATCH --job-name overlap
#SBATCH --account=dpicls
#SBATCH -o /gscratch/ljin1/data/twitter/overlap.%j.out
#SBATCH -e /gscratch/ljin1/data/twitter/overlap.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ljin1@uwyo.edu
#SBATCH --time=6-23:59:59

#################### teton-knl #######################
#SBATCH --partition=teton-knl
#SBATCH --nodes=12
#SBATCH --ntasks-per-node=288
#SBATCH --mem-per-cpu=1300

# #################### teton-hugemem #######################
# #SBATCH --partition=teton-hugemem
# #SBATCH --nodes=10
# #SBATCH --ntasks-per-node=32
# #SBATCH --mem-per-cpu=32000
#
# #################### teton-gpu #######################
# #SBATCH --partition=teton-gpu
# #SBATCH --nodes=8
# #SBATCH --ntasks-per-node=32
# #SBATCH --mem-per-cpu=16000
#
# #################### teton #######################
# #SBATCH --partition=teton
# #SBATCH --nodes=8
# #SBATCH --ntasks-per-node=32
# #SBATCH --mem-per-cpu=4000
#
# #################### moran #######################
# #SBATCH --partition=moran
# #SBATCH --nodes=16
# #SBATCH --ntasks-per-node=16
# #SBATCH --mem-per-cpu=8000
#
# #################### moran-bigmem-gpu #######################
# #SBATCH --partition=moran-bigmem-gpu
# #SBATCH --nodes=2
# #SBATCH --ntasks-per-node=16
# #SBATCH --mem-per-cpu=32000
#
# #################### moran-hugemem #######################
# #SBATCH --partition=moran-hugemem
# #SBATCH --nodes=2
# #SBATCH --ntasks-per-node=8
# #SBATCH --mem-per-cpu=128000


# Start the job
module use ~/.modulefiles
if [[ $OSTYPE == darwin* ]] # MacBook Pro @ libaooutrage (macOS)
then
    DATA_DIR=/Users/libao/Document/data/twitter/result
elif [[ $OSTYPE == linux-gnu ]] # Teton @ UWyo ARCC (Linux)
then
    DATA_DIR=/gscratch/ljin1/data/twitter/result
    module load petsc/3.10.3-gnu
fi
# Set up output directory
if [ ! -d $DATA_DIR ]; then
    mkdir -p $DATA_DIR
fi

mpirun -n 2048 python3 ./test_ellipse.py

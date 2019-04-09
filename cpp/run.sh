#!/bin/bash

ARMADILLO_VER=9.300.2

if [[ $OSTYPE == darwin* ]]     # MacBook Pro @ libaooutrage (macOS)
then
    PROJECT_DIR=/Users/libao/Documents/work/projects/research/Tweet-Spatial-Analysis
    source /usr/local/opt/lmod/init/profile
elif [[ $OSTYPE == linux-gnu ]] # Teton @ UWyo ARCC (Linux)
then
    PROJECT_DIR=/home/ljin1/repos/Tweet-Spatial-Analysis
fi

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

# Load modules
module --version
module purge
module use ~/.modulefiles
module load armadillo/$ARMADILLO_VER
module list

# Run
cd $PROJECT_DIR/cpp
bin/build_distance_matrix -o $DATA_DIR/distance_matrix.csv -i $DATA_DIR/adjacency_matrix.csv

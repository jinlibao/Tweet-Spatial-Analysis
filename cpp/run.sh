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

# Load modules
module --version
module purge
module use ~/.modulefiles
module load armadillo/$ARMADILLO_VER
module list

# Run
cd $PROJECT_DIR/cpp
bin/build_distance_matrix
bin/example

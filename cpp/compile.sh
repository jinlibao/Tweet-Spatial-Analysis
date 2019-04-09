#!/bin/bash

ARMADILLO_VER=9.300.2

if [[ $OSTYPE == darwin* ]]     # MacBook Pro @ libaooutrage (macOS)
then
    PROJECT_DIR=/Users/libao/Documents/work/projects/research/Tweet-Spatial-Analysis
    export MKL_INCLUDE=/opt/intel/mkl/include
    export MKL_LIB=/opt/intel/mkl/lib
    export ARMADILLO_INCLUDE=/Users/libao/.armadillo/armadillo-$ARMADILLO_VER-mkl/include
    export ARMADILLO_LIB=/Users/libao/.armadillo/armadillo-$ARMADILLO_VER-mkl/lib
    source /usr/local/opt/lmod/init/profile
elif [[ $OSTYPE == linux-gnu ]] # Teton @ UWyo ARCC (Linux)
then
    PROJECT_DIR=/home/ljin1/repos/Tweet-Spatial-Analysis
fi

# Set up directory for output data
if [ ! -d $PROJECT_DIR/cpp/bin ]; then
    mkdir -p $PROJECT_DIR/cpp/bin
fi

# Load modules
module --version
module purge
module use ~/.modulefiles
module load armadillo/$ARMADILLO_VER
module list

# Compile
cd $PROJECT_DIR/cpp/src
make distclean -s
make all -s
make clean -s
echo "All compiled!"
cd ..
echo "executable file(s) under './bin/'"
ls ./bin/

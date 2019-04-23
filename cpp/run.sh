#!/bin/bash
#
# run.sh
# Copyright (C) 2019 Libao Jin <jinlibao@outlook.com>
#
# Distributed under terms of the MIT license.

ARMADILLO_VER=9.300.2

if [[ $OSTYPE == darwin* ]]     # MacBook Pro @ libaooutrage (macOS)
then
    source /usr/local/opt/lmod/init/profile
    module purge -q
    module use ~/.modulefiles
    module load armadillo/$ARMADILLO_VER-mkl
    PROJECT_DIR=/Users/libao/Documents/work/projects/research/Tweet-Spatial-Analysis
    DATA_DIR=/Users/libao/Documents/data/twitter/csv
    export MKL_INCLUDE=/opt/intel/mkl/include
    export MKL_LIB=/opt/intel/mkl/lib
    export ARMADILLO_ROOT=/Users/libao/.armadillo/armadillo-$ARMADILLO_VER-mkl
    export ARMADILLO_INCLUDE=$ARMADILLO_ROOT/include
    export ARMADILLO_LIB=$ARMADILLO_ROOT/lib
elif [[ $OSTYPE == linux-gnu ]] # Teton @ UWyo ARCC (Linux)
then
    module purge -q
    module use ~/.modulefiles
    module load swset/2018.05 gcc/7.3.0 openmpi/3.1.0 intel-mkl/2018.2.199 armadillo/$ARMADILLO_VER-mkl
    PROJECT_DIR=/home/ljin1/repos/Tweet-Spatial-Analysis
    DATA_DIR=/gscratch/ljin1/data/twitter/csv
    export INTEL_ROOT=/apps/u/gcc/7.3.0/intel-mkl/2018.2.199-isxsqpg
    export MKL_ROOT=$INTEL_ROOT/compilers_and_libraries_2018.2.199/linux
    export MKL_INCLUDE=$MKL_ROOT/mkl/include
    export MKL_LIB=$MKL_ROOT/mkl/lib/intel64_lin
    export LD_LIBRARY_PATH=$MKL_ROOT/tbb/lib/intel64_lin/gcc4.7:$MKL_ROOT/compiler/lib/intel64_lin:$MKL_ROOT/mkl/lib/intel64_lin:$INTEL_ROOT/lib:$LD_LIBRARY_PATH
    export LIBRARY_PATH=$MKL_ROOT/tbb/lib/intel64_lin/gcc4.7:$MKL_ROOT/compiler/lib/intel64_lin:$MKL_ROOT/mkl/lib/intel64_lin:$INTEL_ROOT/lib:$LIBRARY_PATH
    export ARMADILLO_ROOT=/home/ljin1/.armadillo/armadillo-$ARMADILLO_VER-mkl
    export ARMADILLO_INCLUDE=$ARMADILLO_ROOT/include
    export ARMADILLO_LIB=$ARMADILLO_ROOT/lib64
fi

# Set up output directory
if [ ! -d $DATA_DIR ]; then
    mkdir -p $DATA_DIR
fi

# Run
cd $PROJECT_DIR/cpp

ROWS=5790
ELLIPSECSV[0]=$PROJECT_DIR/data/tweets_mean_all_filtered.csv
ADJ_MATRIX[0]=$DATA_DIR/tweets_mean_all_adjacency_matrix.csv
DIS_MATRIX[0]=$DATA_DIR/tweets_mean_all_distance_matrix.csv
ELLIPSECSV[1]=$PROJECT_DIR/data/tweets_median_working_filtered.csv
ADJ_MATRIX[1]=$DATA_DIR/tweets_median_working_adjacency_matrix.csv
DIS_MATRIX[1]=$DATA_DIR/tweets_median_working_distance_matrix.csv
ELLIPSECSV[2]=$PROJECT_DIR/data/tweets_median_non_working_filtered.csv
ADJ_MATRIX[2]=$DATA_DIR/tweets_median_non_working_adjacency_matrix.csv
DIS_MATRIX[2]=$DATA_DIR/tweets_median_non_working_distance_matrix.csv

for (( i=0; i < 3; ++i ))
do
    bin/main -r $ROWS -e ${ELLIPSECSV[i]} -i ${ADJ_MATRIX[i]} -o ${DIS_MATRIX[i]}
done

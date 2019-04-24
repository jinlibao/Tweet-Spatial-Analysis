#!/bin/bash

#SBATCH --job-name overlap
#SBATCH --account=dpicls
#SBATCH --output=/gscratch/ljin1/data/twitter/log/overlap.%j.out
#SBATCH --error=/gscratch/ljin1/data/twitter/log/overlap.%j.err
#SBATCH --mail-type=END
#SBATCH --mail-user=ljin1@uwyo.edu
#SBATCH --time=6-23:59:59

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
    module load arcc/0.1 slurm/18.08 swset/2018.05 gcc/7.3.0 openmpi/3.1.0 intel-mkl/2018.2.199 armadillo/$ARMADILLO_VER-mkl
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

if [ -z $ELLIPSECSV_FILE ]; then
    ELLIPSECSV_FILE=$PROJECT_DIR/data/tweets_median_working_filtered.csv
    ADJ_MATRIX_FILE=$DATA_DIR/tweets_median_working_adjacency_matrix.csv
    ADJ_ORDER_FILE=$DATA_DIR/tweets_median_working_adjacency_matrix_ordered.csv
    DIS_MATRIX_FILE=$DATA_DIR/tweets_median_working_distance_matrix.csv
fi

if [ -z $ROWS ]; then
    ROWS=500
fi

bin/main -r $ROWS -e $ELLIPSECSV_FILE -a $ADJ_MATRIX_FILE -o $ADJ_ORDER_FILE -d $DIS_MATRIX_FILE

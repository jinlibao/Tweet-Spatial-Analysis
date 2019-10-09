#! /bin/bash
#
# install_armadillo.sh
# Copyright (C) 2019 Libao Jin <jinlibao@outlook.com>
#
# Distributed under terms of the MIT license.

ARMADILLO_VER=9.800.1
ROOT=~/.armadillo
LOC=armadillo-$ARMADILLO_VER
DIR=$ROOT/$LOC

if [[ $OSTYPE == darwin* ]] # MacBook Pro @ libaooutrage (macOS)
then
    export LD_LIBRARY_PATH=/opt/intel/mkl/lib:$LD_LIBRARY_PATH
    export MKL_ROOT=/opt/intel/mkl
    export MKL_INCLUDE=$MKL_ROOT/include
    export MKL_LIBRARY=$MKL_ROOT/lib
    source $MKL_ROOT/bin/mklvars.sh intel64
    source /opt/intel/bin/compilervars.sh intel64
    export CMAKE_INCLUDE_PATH=$MKL_INCLUDE:$CMAKE_INCLUDE_PATH
    export CMAKE_LIBRARY_PATH=$MKL_LIBRARY:$CMAKE_LIBRARY_PATH
elif [[ $OSTYPE == linux-gnu ]] # Teton @ UWyo ARCC (Linux)
then
    module purge -q
    module use ~/.modulefiles
    module load swset/2018.05 gcc/7.3.0 openmpi/3.1.0 intel-mkl/2018.2.199
    PROJECT_DIR=/home/ljin1/repos/Tweet-Spatial-Analysis
    export INTEL_ROOT=/apps/u/gcc/7.3.0/intel-mkl/2018.2.199-isxsqpg
    export MKL_ROOT=$INTEL_ROOT/compilers_and_libraries_2018.2.199/linux
    export MKL_INCLUDE=$MKL_ROOT/mkl/include
    export MKL_LIB=$MKL_ROOT/mkl/lib/intel64_lin
    export LD_LIBRARY_PATH=$MKL_ROOT/tbb/lib/intel64_lin/gcc4.7:$MKL_ROOT/compiler/lib/intel64_lin:$MKL_ROOT/mkl/lib/intel64_lin:$INTEL_ROOT/lib:$LD_LIBRARY_PATH
    export LIBRARY_PATH=$MKL_ROOT/tbb/lib/intel64_lin/gcc4.7:$MKL_ROOT/compiler/lib/intel64_lin:$MKL_ROOT/mkl/lib/intel64_lin:$INTEL_ROOT/lib:$LIBRARY_PATH
    export CMAKE_INCLUDE_PATH=$MKL_INCLUDE:$CMAKE_INCLUDE_PATH
    export CMAKE_LIBRARY_PATH=$LD_LIBRARY_PATH:$LIBRARY_PATH:$MKL_LIBRARY:$CMAKE_LIBRARY_PATH
    source $INTEL_ROOT/bin/compilervars.sh intel64
    source $MKL_ROOT/bin/compilervars.sh intel64
    source $MKL_ROOT/mkl/bin/mklvars.sh intel64
fi

if [ ! -d $ROOT ]; then
    mkdir -p $ROOT
fi

cd $ROOT

if [ ! -f $LOC.tar.xz ]; then
    wget http://sourceforge.net/projects/arma/files/armadillo-$ARMADILLO_VER.tar.xz
fi

if [ ! -d $LOC ]; then
    tar xf $LOC.tar.xz
fi

cd $LOC

if [ -d $DIR-mkl ]; then
    rm -rf $DIR-mkl
fi

CC=gcc CXX=g++ cmake . -DCMAKE_INSTALL_PREFIX:PATH=$DIR-mkl
make
make install

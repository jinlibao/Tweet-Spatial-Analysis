#! /bin/sh
#
# install_armadillo.sh
# Copyright (C) 2019 Libao Jin <jinlibao@outlook.com>
#
# Distributed under terms of the MIT license.
#


ARMADILLO_VER=9.300.2
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
    export LD_LIBRARY_PATH=/opt/intel/mkl/lib:$LD_LIBRARY_PATH
    export MKL_ROOT=/opt/intel/mkl
    export MKL_INCLUDE=$MKL_ROOT/include
    export MKL_LIBRARY=$MKL_ROOT/lib
    source $MKL_ROOT/bin/mklvars.sh intel64
    source /opt/intel/bin/compilervars.sh intel64
    export CMAKE_INCLUDE_PATH=$MKL_INCLUDE:$CMAKE_INCLUDE_PATH
    export CMAKE_LIBRARY_PATH=$MKL_LIBRARY:$CMAKE_LIBRARY_PATH
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

cmake . -DCMAKE_INSTALL_PREFIX:PATH=$DIR-mkl
make
make install


#! /bin/bash
#
# submit_task.sh
# Copyright (C) 2019 Libao Jin <jinlibao@outlook.com>
#
# Distributed under terms of the MIT license.
#

MORAN_PARTITION=moran
MORAN_MAX_NODES=128
MORAN_NODES=16
MORAN_NTASKS_PER_NODE=16
MORAN_MEM=128000

MORAN_GPU_PARTITION=moran-bigmem-gpu
MORAN_GPU_MAX_NODES=2
MORAN_GPU_NODES=2
MORAN_GPU_NTASKS_PER_NODE=16
MORAN_GPU_MEM=512000

MORAN_HUGEMEM_PARTITION=moran-hugemem
MORAN_HUGEMEM_MAX_NODES=2
MORAN_HUGEMEM_NODES=2
MORAN_HUGEMEM_NTASKS_PER_NODE=8
MORAN_HUGEMEM_MEM=1024000

TETON_PARTITION=teton
TETON_MAX_NODES=128
#TETON_NODES=8
#TETON_NTASKS_PER_NODE=32
TETON_NODES=1
TETON_NTASKS_PER_NODE=1
TETON_MEM=128000

TETON_GPU_PARTITION=teton-gpu
TETON_GPU_MAX_NODES=8
TETON_GPU_NODES=8
TETON_GPU_NTASKS_PER_NODE=32
TETON_GPU_MEM=512000

TETON_HUGEMEM_PARTITION=teton-hugemem
TETON_HUGEMEM_MAX_NODES=10
TETON_HUGEMEM_NODES=10
TETON_HUGEMEM_NTASKS_PER_NODE=32
TETON_HUGEMEM_MEM=1024000

TETON_KNL_PARTITION=teton-knl
TETON_KNL_MAX_NODES=12
TETON_KNL_NODES=11
TETON_KNL_NTASKS_PER_NODE=288
TETON_KNL_MEM=384000

PARTITION[0]=$MORAN_PARTITION
NODES[0]=$MORAN_NODES
NTASKS_PER_NODE[0]=$MORAN_NTASKS_PER_NODE
MEM[0]=$MORAN_MEM

PARTITION[1]=$MORAN_GPU_PARTITION
NODES[1]=$MORAN_GPU_NODES
NTASKS_PER_NODE[1]=$MORAN_GPU_NTASKS_PER_NODE
MEM[1]=$MORAN_GPU_MEM

PARTITION[2]=$MORAN_HUGEMEM_PARTITION
NODES[2]=$MORAN_HUGEMEM_NODES
NTASKS_PER_NODE[2]=$MORAN_HUGEMEM_NTASKS_PER_NODE
MEM[2]=$MORAN_HUGEMEM_MEM

PARTITION[3]=$TETON_PARTITION
NODES[3]=$TETON_NODES
NTASKS_PER_NODE[3]=$TETON_NTASKS_PER_NODE
MEM[3]=$TETON_MEM

PARTITION[4]=$TETON_GPU_PARTITION
NODES[4]=$TETON_GPU_NODES
NTASKS_PER_NODE[4]=$TETON_GPU_NTASKS_PER_NODE
MEM[4]=$TETON_GPU_MEM

PARTITION[5]=$TETON_HUGEMEM_PARTITION
NODES[5]=$TETON_HUGEMEM_NODES
NTASKS_PER_NODE[5]=$TETON_HUGEMEM_NTASKS_PER_NODE
MEM[5]=$TETON_HUGEMEM_MEM

PARTITION[6]=$TETON_KNL_PARTITION
NODES[6]=$TETON_KNL_NODES
NTASKS_PER_NODE[6]=$TETON_KNL_NTASKS_PER_NODE
MEM[6]=$TETON_KNL_MEM

if [[ $OSTYPE == darwin* ]]     # MacBook Pro @ libaooutrage (macOS)
then
    PROJECT_DIR=/Users/libao/Documents/work/projects/research/Tweet-Spatial-Analysis
    DATA_DIR=/Users/libao/Documents/data/twitter/csv
elif [[ $OSTYPE == linux-gnu ]] # Teton @ UWyo ARCC (Linux)
then
    PROJECT_DIR=/home/ljin1/repos/Tweet-Spatial-Analysis
    DATA_DIR=/gscratch/ljin1/data/twitter/csv
fi

# Set up output directory
if [ ! -d $DATA_DIR ]; then
    mkdir -p $DATA_DIR
fi

k=3
ROW[0]=2000
ROW[1]=4000
ROW[2]=8000
ROW[3]=16000
ROW[4]=57909

ELLIPSECSV[0]=$PROJECT_DIR/data/tweets_mean_all_filtered.csv
ELLIPSECSV[1]=$PROJECT_DIR/data/tweets_median_working_filtered.csv
ELLIPSECSV[2]=$PROJECT_DIR/data/tweets_median_non_working_filtered.csv
ADJ_MATRIX[0]=$DATA_DIR/tweets_mean_all_adjacency_matrix.csv
ADJ_MATRIX[1]=$DATA_DIR/tweets_median_working_adjacency_matrix.csv
ADJ_MATRIX[2]=$DATA_DIR/tweets_median_non_working_adjacency_matrix.csv
DIS_MATRIX[0]=$DATA_DIR/tweets_mean_all_distance_matrix.csv
DIS_MATRIX[1]=$DATA_DIR/tweets_median_working_distance_matrix.csv
DIS_MATRIX[2]=$DATA_DIR/tweets_median_non_working_distance_matrix.csv
 ADJ_ORDER[0]=$DATA_DIR/tweets_mean_all_adjacency_matrix_ordered.csv
 ADJ_ORDER[1]=$DATA_DIR/tweets_median_working_adjacency_matrix_ordered.csv
 ADJ_ORDER[2]=$DATA_DIR/tweets_median_non_working_adjacency_matrix_ordered.csv

for (( j=4; j<5; ++j ))
do
    export ROWS=${ROW[j]}
    for (( i=0; i < 3; ++i ))
    do
        export ELLIPSECSV_FILE=${ELLIPSECSV[i]}
        export ADJ_MATRIX_FILE=${ADJ_MATRIX[i]}
        export ADJ_ORDER_FILE=${ADJ_ORDER[i]}
        export DIS_MATRIX_FILE=${DIS_MATRIX[i]}
        export NCPU=`echo ${NODES[k]} \* ${NTASKS_PER_NODE[k]}|bc`
        sbatch --partition=${PARTITION[k]} --nodes=${NODES[k]} --ntasks-per-node=${NTASKS_PER_NODE[k]} --mem=${MEM[k]} task.sh
    done
done

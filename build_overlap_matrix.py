#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Libao Jin <jinlibao@outlook.com>
#
# Distributed under terms of the MIT license.

"""
Build overlap matrix
"""

import math
import timeit
import numpy as np
import pandas as pd
from utils import analysis_utilities as au
from utils import configuration_utilities as cu
from utils import file_utilities as fu
from utils import tweet_data_utilities as tdu
from mpi4py import MPI

def build_overlap_matrix(df, filename='overlap_matrix.csv'):
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    if rank == 0:
        # rows, cols = df.shape[0], df.shape[0]
        rows, cols = 20, 20
        tweet_data_working_overlap = np.zeros((rows, rows))

        for i in range(rows):
            for j in range(i, rows):
                if i == j:
                    tweet_data_working_overlap[i][j] = 1
                    continue
                elif df['a'][i] == 0 or df['b'][i] == 0 or df['a'][j] == 0 or df['b'][j] == 0:
                    tweet_data_working_overlap[i][j] = -1
                    continue
                tweet_data_working_overlap[i][j] = au.are_two_ellipses_overlapping(
                    df['x'][i], df['y'][i], df['a'][i], df['b'][i], df['angle'][i],
                    df['x'][j], df['y'][j], df['a'][j], df['b'][j], df['angle'][j],
                )
        for i in range(rows):
            for j in range(0, i):
                tweet_data_working_overlap[i][j] = tweet_data_working_overlap[j][i]
        tweet_data_working_overlap = pd.DataFrame(data=tweet_data_working_overlap.astype(int), columns=df['id'][0:rows], index=df['id'][0:rows])
        # print(tweet_data_working_overlap)
        tweet_data_working_overlap.to_csv(filename, sep=',', header=True, index=True)

def build_overlap_matrix_parallel(df, filename='overlap_matrix.csv'):
    start_time = timeit.default_timer()
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    # rows, cols = df.shape
    rows, cols = 20, 10
    rows_per_cpu = math.ceil(rows / size)
    k = rank

    if k != 0:
        start = k * rows_per_cpu
        if k == size - 1:
            # rows_local = (rows - 1) % rows_per_cpu + 1
            rows_local = rows - k * rows_per_cpu
            if rows_local < 0:
                rows_local = 0
        else:
            rows_local = rows_per_cpu
        tweet_data_working_overlap = np.zeros((rows_local, rows))
        for i in range(0, rows_local):
            for j in range(rows):
                if df['a'][start + i] == 0 or df['b'][start + i] == 0 or df['a'][j] == 0 or df['b'][j] == 0:
                    tweet_data_working_overlap[i][j] = -1
                    continue
                tweet_data_working_overlap[i][j] = au.are_two_ellipses_overlapping(
                    df['x'][start + i], df['y'][start + i], df['a'][start + i], df['b'][start + i], df['angle'][start + i],
                    df['x'][j], df['y'][j], df['a'][j], df['b'][j], df['angle'][j],
                )
        # print(tweet_data_working_overlap, rank)
        # print('My rank is {:d} of {:d} processors, number of rows: {:d}, and rows_per_cpu: {:d}'.format(rank, size, rows, rows_local))
        comm.send(tweet_data_working_overlap, dest=0)

    if k == 0:
        data_dict = {}
        start = k * rows_per_cpu
        if k == size - 1:
            rows_local = (rows - 1) % rows_per_cpu + 1
        else:
            rows_local = rows_per_cpu
        tweet_data_working_overlap = np.zeros((rows_local, rows))
        for i in range(0, rows_local):
            for j in range(rows):
                if df['a'][start + i] == 0 or df['b'][start + i] == 0 or df['a'][j] == 0 or df['b'][j] == 0:
                    tweet_data_working_overlap[i][j] = -1
                    continue
                tweet_data_working_overlap[i][j] = au.are_two_ellipses_overlapping(
                    df['x'][start + i], df['y'][start + i], df['a'][start + i], df['b'][start + i], df['angle'][start + i],
                    df['x'][j], df['y'][j], df['a'][j], df['b'][j], df['angle'][j],
                )
        # print(tweet_data_working_overlap, rank)
        # print('My rank is {:d} of {:d} processors, number of rows: {:d}, and rows_per_cpu: {:d}'.format(rank, size, rows, rows_local))

        for i in range(1, size):
            data_dict[str(i)] = comm.recv(source=i)
        for i in range(1, len(data_dict) + 1):
            # print(data_dict[str(i)], i)
            # tweet_data_working_overlap.extend(data_dict[str(i)])
            tweet_data_working_overlap = np.append(tweet_data_working_overlap , data_dict[str(i)], axis = 0)
        tweet_data_working_overlap = pd.DataFrame(data=tweet_data_working_overlap.astype(int), columns=df['id'][0:rows], index=df['id'][0:rows])
        # print(tweet_data_working_overlap)
        tweet_data_working_overlap.to_csv(filename, sep=',', header=True, index=True)
        end_time = timeit.default_timer()
        print(end_time - start_time)

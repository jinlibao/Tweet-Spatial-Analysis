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
import datetime
import platform
import numpy as np
import pandas as pd
from argparse import ArgumentParser
from utils import analysis_utilities as au
from utils import configuration_utilities as cu
from utils import file_utilities as fu
from utils import tweet_data_utilities as tdu
from data_preprocessing import TweetDataPreProcessing
from mpi4py import MPI

def build_overlap_matrix(df, rows, filename='./overlap_matrix.csv'):
    start_time = timeit.default_timer()
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    if rank == 0:
        tweet_data_working_overlap = np.zeros((rows, rows))
        for i in range(rows):
            for j in range(rows):
                if i == j:
                    tweet_data_working_overlap[i, j] = 1
                elif j < i:
                    if df['a'][i] == 0 or df['b'][i] == 0 or df['a'][j] == 0 or df['b'][j] == 0:
                        tweet_data_working_overlap[i, j] = -1
                    else:
                        tweet_data_working_overlap[i, j] = au.are_two_ellipses_overlapping(
                            df['x'][i], df['y'][i], df['a'][i], df['b'][i], df['angle'][i],
                            df['x'][j], df['y'][j], df['a'][j], df['b'][j], df['angle'][j],
                        )
            if (i + 1) % 10 == 0:
                print('CPU {:04d}: {:6.2f}% accomplished'.format(rank,  (i + 1) / rows * 100))

        for i in range(rows):
            tweet_data_working_overlap[0:i, i] = tweet_data_working_overlap[i, 0:i]
        print('CPU {:04d}: Assigned transposed tril to triu (Time: {:.2f} seconds)'.format(rank, timeit.default_timer() - start_time))

        tweet_data_working_overlap = pd.DataFrame(data=tweet_data_working_overlap.astype(int), columns=df['id'][0:rows], index=df['id'][0:rows])
        print('CPU {:04d}: Converted numpy array to pandas DataFrame (Time: {:.2f} seconds)'.format(rank, timeit.default_timer() - start_time))
        # print(tweet_data_working_overlap)
        # print('CPU {:04d}: Now saving pandas DataFrame to {:s} (Time: {:.2f} seconds)'.format(rank, filename, timeit.default_timer() - start_time))
        # tweet_data_working_overlap.to_csv(filename, sep=',', header=True, index=True)
        bdm.build_distance_matrix(tweet_data_working_overlap, filename, start_time)

        elapsed_time = timeit.default_timer() - start_time
        hour = math.floor(elapsed_time / 3600)
        minute = math.floor((elapsed_time - hour * 3600) / 60)
        second = elapsed_time - 3600 * hour - 60 * minute
        print('Time elapsed: {:d} hours {:d} minutes {:.2f} seconds'.format(hour, minute, second))
        return tweet_data_working_overlap

def build_overlap_matrix_parallel(df, rows, filename='./overlap_matrix_parallel.csv'):
    start_time = timeit.default_timer()
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    status = MPI.Status()
    row_start, rows_local = find_range_1D(rows, size, rank)
    # print('Total CPUs: {:4d}, Rank {:4d}, rows: {:5d}, rows_per_cpu: {:5d}, row_start: {:5d}'.format(size, rank, rows, rows_local, row_start))

    tweet_data_working_overlap_local = np.zeros((rows_local, rows))
    for i in range(row_start, row_start + rows_local):
        for j in range(rows):
            if i == j:
                tweet_data_working_overlap_local[i - row_start, j] = 1
            elif j < i:
                if df['a'][i] == 0 or df['b'][i] == 0 or df['a'][j] == 0 or df['b'][j] == 0:
                    tweet_data_working_overlap_local[i - row_start, j] = -1
                else:
                    tweet_data_working_overlap_local[i - row_start, j] = au.are_two_ellipses_overlapping(
                        df['x'][i], df['y'][i], df['a'][i], df['b'][i], df['angle'][i],
                        df['x'][j], df['y'][j], df['a'][j], df['b'][j], df['angle'][j],
                    )

        if (i - row_start + 1) % 10 == 0:
            print('CPU {:04d}: {:6.2f}% accomplished'.format(rank,  (i - row_start + 1) / rows_local * 100))
    # print(tweet_data_working_overlap, rank)

    if rank != 0:
        comm.send(tweet_data_working_overlap_local, dest=0)

    if rank == 0:
        tweet_data_working_overlap = np.zeros((rows, rows))
        tweet_data_working_overlap[0:rows_local] = tweet_data_working_overlap_local
        for i in range(row_start, row_start + rows_local):
            tweet_data_working_overlap[0:i, i] = tweet_data_working_overlap[i, 0:i]
        number_of_sources = 1
        while number_of_sources < size:
            data = comm.recv(source=MPI.ANY_SOURCE, status=status)
            k = status.Get_source()
            print('CPU {:04d}: Received data from CPU {:04d} (Time: {:.2f} seconds)'.format(rank,  k, timeit.default_timer() - start_time))
            row_start, rows_local = find_range_1D(rows, size, k)
            tweet_data_working_overlap[row_start:row_start + rows_local] = data + tweet_data_working_overlap[row_start:row_start + rows_local]
            for i in range(row_start, row_start + rows_local):
                tweet_data_working_overlap[0:i, i] = tweet_data_working_overlap[i, 0:i]
            number_of_sources += 1
        print('CPU {:04d}: Merged received data to global numpy array (Time: {:.2f} seconds)'.format(rank, timeit.default_timer() - start_time))

        tweet_data_working_overlap = pd.DataFrame(data=tweet_data_working_overlap.astype(int), columns=df['id'][0:rows], index=df['id'][0:rows])
        print('CPU {:04d}: Converted numpy array to pandas DataFrame (Time: {:.2f} seconds)'.format(rank, timeit.default_timer() - start_time))
        # print(tweet_data_working_overlap)
        # print('CPU {:04d}: Now saving pandas DataFrame to {:s} (Time: {:.2f} seconds)'.format(rank, filename, timeit.default_timer() - start_time))
        # tweet_data_working_overlap.to_csv(filename, sep=',', header=True, index=True)
        bdm.build_distance_matrix(tweet_data_working_overlap, filename, start_time)

        elapsed_time = timeit.default_timer() - start_time
        hour = math.floor(elapsed_time / 3600)
        minute = math.floor((elapsed_time - hour * 3600) / 60)
        second = elapsed_time - 3600 * hour - 60 * minute
        print('Time elapsed: {:d} hours {:d} minutes {:.2f} seconds'.format(hour, minute, second))
        return tweet_data_working_overlap

def bfs_array(raw_graph):
    '''
    Breadth First Search on array
    Find components of a graph
    '''
    visited = [False for i in range(len(raw_graph))]
    components = []
    for i in range(len(raw_graph)):
        if not visited[i]:
            visited[i] = True
            component = [i]
            queue = []
            for j in range(len(raw_graph[i])):
                if not visited[j] and raw_graph[i, j] > 0:
                    visited[j] = True
                    queue.append(j)
                    component.append(j)
            while len(queue) > 0:
                j = queue.pop(0)
                for k in range(len(raw_graph[j])):
                    if not visited[k] and raw_graph[j, k] > 0:
                        visited[k] = True
                        component.append(k)
                        queue.append(k)
            components.append(component)
    return components

def bfs_df(raw_graph):
    '''
    Breadth First Search on pandas.DataFrame
    Find components of a graph
    '''
    n, _ = raw_graph.shape
    visited = [False for i in range(n)]
    components = []
    for i in range(n):
        if not visited[i]:
            visited[i] = True
            component = [i]
            queue = []
            for j in range(n):
                if not visited[j] and raw_graph.iloc[i, j] > 0:
                    visited[j] = True
                    queue.append(j)
                    component.append(j)
            while len(queue) > 0:
                j = queue.pop(0)
                for k in range(n):
                    if not visited[k] and raw_graph.iloc[j, k] > 0:
                        visited[k] = True
                        component.append(k)
                        queue.append(k)
            components.append(component)
    return components

def APD_recursive(A):
    '''
    Seidel's algorithm
    Find distance matrix of a connected graph
    '''
    n, _ = A.shape
    if all(A[i,j] for i in range(n) for j in range(n) if i != j): return A
    Z = A.dot(A)
    B = np.array([[1 if i != j and ( A[i,j] == 1 or Z[i,j] > 0 ) else 0 for j in range(n)] for i in range(n)])
    T = APD_recursive(B)
    X = T.dot(A)
    degree = [sum( A[i,j] for j in range(n) ) for i in range(n)]
    D = np.array([[2 * T[i,j] if X[i,j] >= T[i,j] * degree[j] else 2 * T[i,j] - 1 for j in range(n)] for i in range(n) ] )
    return D

def APD(A):
    '''
    Seidel's algorithm
    Find distance matrix of a connected graph
    '''
    n, _ = A.shape
    B = A
    while not all(B[i,j] for i in range(n) for j in range(n) if i != j):
        Z = B.dot(B)
        B = np.array([[1 if i != j and ( B[i,j] == 1 or Z[i,j] > 0 ) else 0 for j in range(n)] for i in range(n)])
    T = B
    X = T.dot(A)
    degree = [sum( A[i,j] for j in range(n) ) for i in range(n)]
    D = np.array([[2 * T[i,j] if X[i,j] >= T[i,j] * degree[j] else 2 * T[i,j] - 1 for j in range(n)] for i in range(n) ] )
    return D

def find_components(df, filename, start_time):
    print('CPU {:04d}: Now cleaning data remove -1 entries (Time: {:.2f} seconds)'.format(0, timeit.default_timer() - start_time))
    # Data cleaning
    col = [idx for idx in df.iloc[0].index if df.iloc[0][idx] == -1]
    idx = [int(idx) for idx in col]
    df = df.drop(index=idx, columns=col)

    print('CPU {:04d}: Starting BFS... (Time: {:.2f} seconds)'.format(0, timeit.default_timer() - start_time))
    components = bfs_df(df)
    components.sort(key=len, reverse=True)

    idx = df.index[components[0]]
    col = df.columns[components[0]]
    idx_list = [idx]
    col_list = [col]
    for i in range(1, len(components)):
        idx = idx.append(df.index[components[i]])
        idx_list.append(df.index[components[i]])
        col = col.append(df.columns[components[i]])
        col_list.append(df.columns[components[i]])
        print(components[i], idx.shape, col.shape, type(idx), type(col))

    print('CPU {:04d}: Reconstructing adjacency pandas DataFrame (Time: {:.2f} seconds)'.format(0, timeit.default_timer() - start_time))
    df_orig = np.zeros(df.shape)
    df_orig = pd.DataFrame(data=df_orig.astype(int), columns=col, index=idx)
    for i in range(len(idx_list)):
        df_orig.loc[idx_list[i], col_list[i]] = df.loc[idx_list[i], col_list[i]]
    filename_adj = filename.replace('.csv', '_adjacency.csv')
    print('CPU {:04d}: Now saving pandas DataFrame to {:s} (Time: {:.2f} seconds)'.format(0, filename_adj, timeit.default_timer() - start_time))
    df_orig.to_csv(filename_adj, sep=',', header=True, index=True)
    return (df, components, idx, col, idx_list, col_list)

def build_distance_matrix(df, components, idx, col, idx_list, col_list, filename, start_time):
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    status = MPI.Status()

    df_dist_list = []
    for i in range(1, len(components)):
        if rank == i % size:
            print('CPU {:04d}: Starting APD... (Time: {:.2f} seconds)'.format(rank, timeit.default_timer() - start_time))
            D = APD_recursive(np.array(df.loc[idx_list[i],col_list[i]]))
            comm.send((D, i), dest=0)

    if rank == 0:
        print('CPU {:04d}: Starting APD... (Time: {:.2f} seconds)'.format(rank, timeit.default_timer() - start_time))
        D = APD_recursive(np.array(df.loc[idx_list[rank],col_list[rank]]))
        df_dist_list.append(pd.DataFrame(data=D.astype(int), columns=col_list[rank], index=idx_list[rank]))
        number_of_sources = 1

        while number_of_sources < len(components):
            D, i = comm.recv(source=MPI.ANY_SOURCE, status=status)
            k = status.Get_source()
            print('CPU {:04d}: Received data from CPU {:04d} (Time: {:.2f} seconds)'.format(rank,  k, timeit.default_timer() - start_time))
            df_dist_list.append(pd.DataFrame(data=D.astype(int), columns=col_list[i], index=idx_list[i]))
            number_of_sources += 1

        print('CPU {:04d}: Constructing distance pandas DataFrame (Time: {:.2f} seconds)'.format(0, timeit.default_timer() - start_time))
        df_dist = np.ones(df.shape) * -1
        df_dist = pd.DataFrame(data=df_dist.astype(int), columns=col, index=idx)

        for i in range(len(idx_list)):
            df_dist.loc[idx_list[i], col_list[i]] = df_dist_list[i]

        filename_dis = filename.replace('.csv', '_distance.csv')

        print('CPU {:04d}: Now saving pandas DataFrame to {:s} (Time: {:.2f} seconds)'.format(0, filename_dis, timeit.default_timer() - start_time))
        df_dist.to_csv(filename_dis, sep=',', header=True, index=True)

        elapsed_time = timeit.default_timer() - start_time
        hour = math.floor(elapsed_time / 3600)
        minute = math.floor((elapsed_time - hour * 3600) / 60)
        second = elapsed_time - 3600 * hour - 60 * minute
        print('Time elapsed: {:d} hours {:d} minutes {:.2f} seconds'.format(hour, minute, second))

def build_overlap_matrix_parallel_block(df, rows, filename='./overlap_matrix_parallel_block.csv'):
    start_time = timeit.default_timer()
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    status = MPI.Status()
    cpu_required, row_start, col_start, rows_local, cols_local = find_range_2D(rows, size, rank)
    # print('Total CPUs: {:4d}, Rank {:4d}, rows: {:5d}, rows_per_cpu: {:5d}, row_start: {:5d}, cols_per_cpu: {:5d}, col_start: {:5d}'.format(size, rank, rows, rows_local, row_start, cols_local, col_start))

    if rank < cpu_required:
        tweet_data_working_overlap_local = np.zeros((rows_local, cols_local))
        for i in range(row_start, row_start + rows_local):
            for j in range(col_start, col_start + cols_local):
                if i == j:
                    tweet_data_working_overlap_local[i - row_start, j - col_start] = 0
                elif j < i:
                    if df['a'][i] == 0 or df['b'][i] == 0 or df['a'][j] == 0 or df['b'][j] == 0:
                        tweet_data_working_overlap_local[i - row_start, j - col_start] = -1
                    else:
                        tweet_data_working_overlap_local[i - row_start, j - col_start] = au.are_two_ellipses_overlapping(
                            df['x'][i], df['y'][i], df['a'][i], df['b'][i], df['angle'][i],
                            df['x'][j], df['y'][j], df['a'][j], df['b'][j], df['angle'][j]
                        )
            if row_start == col_start:
                tweet_data_working_overlap_local[0:i - row_start + 1, i - row_start] = tweet_data_working_overlap_local[i - row_start, 0:i - row_start + 1]

            if (i - row_start + 1) % 10 == 0:
                print('CPU {:04d}: {:6.2f}% accomplished'.format(rank,  (i - row_start + 1) / rows_local * 100))

        # print(tweet_data_working_overlap, rank)

        if rank != 0:
            comm.send(tweet_data_working_overlap_local, dest=0)
            # print('CPU {:04d}: Sent data to CPU {:04d}'.format(rank,  0))

        if rank == 0:
            tweet_data_working_overlap = np.zeros((rows, rows))
            tweet_data_working_overlap[0:rows_local, 0:cols_local] = tweet_data_working_overlap_local

            number_of_sources = 1
            while number_of_sources < cpu_required:
                data = comm.recv(source=MPI.ANY_SOURCE, status=status)
                k = status.Get_source()
                print('CPU {:04d}: Received data from CPU {:04d} (Time: {:.2f} seconds)'.format(rank,  k, timeit.default_timer() - start_time))
                cpu_required, row_start, col_start, rows_local, cols_local = find_range_2D(rows, size, k)
                tweet_data_working_overlap[row_start:row_start + rows_local, col_start:col_start + cols_local] = data
                if row_start != col_start:
                    tweet_data_working_overlap[col_start:col_start + cols_local, row_start:row_start + rows_local] = data.T
                number_of_sources += 1
            print('CPU {:04d}: Merged received data to global numpy array (Time: {:.2f} seconds)'.format(rank, timeit.default_timer() - start_time))

            tweet_data_working_overlap = pd.DataFrame(data=tweet_data_working_overlap.astype(int), columns=df['id'][0:rows], index=df['id'][0:rows])
            print('CPU {:04d}: Converted numpy array to pandas DataFrame (Time: {:.2f} seconds)'.format(rank, timeit.default_timer() - start_time))
            # print(tweet_data_working_overlap)
            # print('CPU {:04d}: Now saving pandas DataFrame to {:s} (Time: {:.2f} seconds)'.format(rank, filename, timeit.default_timer() - start_time))
            # tweet_data_working_overlap.to_csv(filename.replace('.csv', '_adjacency.csv'), sep=',', header=True, index=True)
            data = find_components(tweet_data_working_overlap, filename, start_time)

            # if len(components) < size:
                # for i in range(1, len(components)):
                    # comm.send(data, dest=i)

            df, components, idx, col, idx_list, col_list = data
            build_distance_matrix(df, components, idx, col, idx_list, col_list, filename, start_time)

        # if rank != 0:
            # data = comm.recv(source=0, status)
            # df, components, idx, col, idx_list, col_list = data
            # build_distance_matrix(df, components, idx, col, idx_list, col_list, filename, start_time)

def find_range_1D(number, number_of_cpu, rank):
    number_per_cpu = number / number_of_cpu
    number_per_cpu_ceil = math.ceil(number_per_cpu)
    number_per_cpu_floor = math.floor(number_per_cpu)
    threshold = number - number_per_cpu_floor * number_of_cpu

    if rank < threshold:
        start = rank * number_per_cpu_ceil
        number_local = number_per_cpu_ceil
    else:
        start = threshold * number_per_cpu_ceil + (rank - threshold) * number_per_cpu_floor
        number_local = number_per_cpu_floor

    return (start, number_local)

def find_range_2D(rows, cpu_available, rank):
    cols = rows
    m = math.floor(math.sqrt(2 * cpu_available + 1 / 4) - 1 / 2)
    cpu_required = int(m * (m + 1) / 2)
    i = math.floor(math.sqrt(2 * rank + 1 / 4) - 1 / 2)
    j = int(rank - i * (i + 1) / 2)

    row_start, rows_local = find_range_1D(rows, m, i)
    col_start, cols_local = find_range_1D(cols, m, j)

    return (cpu_required, row_start, col_start, rows_local, cols_local)


if __name__ == '__main__':
    file_open = fu.FileOpen("data", "tweet-data.csv")
    tweet_spatial_analysis_config = cu.TweetSpatialAnalysisConfig("conf/tweet_spatial_analysis.ini")
    tdp = TweetDataPreProcessing(file_open, tweet_spatial_analysis_config)
    tdp.read_from_json("data/tweet_mean_all.json",
                       "data/tweets_median_working.json",
                       "data/tweets_median_non_working.json")

    dt = datetime.datetime.now()
    year = dt.year
    month = dt.month
    day= dt.day
    hour = dt.hour
    minute = dt.minute
    second = dt.second

    # algorithm = 'serial'
    # algorithm = 'parallel'
    algorithm = 'parallel_block'

    parser = ArgumentParser()
    parser.add_argument('-l', '--loc', dest='loc', help='Location to save the output file')
    parser.add_argument('-n', '--name', dest='name', help='Name for the output file')
    parser.add_argument('-o', '--outfile', dest='filename', help="Path for saving the output file")
    parser.add_argument('-r', '--rows', dest='rows', help="Path for saving the output file")
    args = parser.parse_args()

    if args.loc:
        loc = args.loc
    else:
        loc = '.'

    if args.name:
        name = args.name
    else:
        if platform.system() == 'Linux':
            node = 'teton'
            # node = 'moran'
            partition = 'hugemem'
            # partition = 'regular'
            # partition = 'knl'
            # partition = 'gpu'
        elif platform.system() == 'Darwin':
            node = 'libaoutrage'
            partition = 'macOS'
        else:
            node = 'intel'
            partition = 'partition'

        name = 'overlap_matrix_{:s}_{:s}_{:s}_{:d}_{:02d}_{:02d}_{:02d}_{:02d}_{:02d}.csv'.format(node, partition, algorithm, year, month, day, hour, minute, second)

    if args.filename:
        filename = args.filename
    else:
        filename = '{:s}/{:s}'.format(loc, name)

    # rows, _ = tdp.tweet_data_working.df.shape
    if args.rows:
        rows = int(args.rows)
    else:
        rows, _ = tdp.tweet_data_working.df.shape
        print(rows)

    if algorithm == 'serial':
        build_overlap_matrix(tdp.tweet_data_working.df, rows, filename)
    elif algorithm == 'parallel':
        build_overlap_matrix_parallel(tdp.tweet_data_working.df, rows, filename)
    else:
        build_overlap_matrix_parallel_block(tdp.tweet_data_working.df, rows, filename)


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
import mpi4py
import numpy as np
import pandas as pd
from argparse import ArgumentParser
from utils import analysis_utilities as au
from utils import configuration_utilities as cu
from utils import file_utilities as fu
from utils import tweet_data_utilities as tdu
from data_preprocessing import TweetDataPreProcessing
from mpi4py import MPI

def build_overlap_matrix(rows, filename='./overlap_matrix.csv'):
    file_open = fu.FileOpen("data", "tweet-data.csv")
    tweet_spatial_analysis_config = cu.TweetSpatialAnalysisConfig("conf/tweet_spatial_analysis.ini")
    tdp = TweetDataPreProcessing(file_open, tweet_spatial_analysis_config)
    tdp.read_from_json("data/tweet_mean_all.json",
                       "data/tweets_median_working.json",
                       "data/tweets_median_non_working.json")
    df = tdp.tweet_data_working.df
    start_time = timeit.default_timer()
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    if rows == 0:
        rows, _ = df.shape
    if rank == 0:
        tweet_data_working_overlap = np.zeros((rows, rows))
        for i in range(rows):
            for j in range(rows):
                if i == j:
                    tweet_data_working_overlap[i, j] = 0
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

        tweet_data_working_overlap = pd.DataFrame(data=tweet_data_working_overlap, columns=df['id'][0:rows], index=df['id'][0:rows])
        print('CPU {:04d}: Converted numpy array to pandas DataFrame (Time: {:.2f} seconds)'.format(rank, timeit.default_timer() - start_time))

        data = find_components(tweet_data_working_overlap, filename, start_time)
        df, components, idx, col, idx_list, col_list = data
        number_of_components = len(components)

        df_dist_list = []
        for i in range(number_of_components):
            D = APD_recursive(np.array(df.loc[idx_list[i],col_list[i]]))
            df_dist_list.append(pd.DataFrame(data=D, columns=col_list[i], index=idx_list[i]))


        print('CPU {:04d}: Constructing distance pandas DataFrame (Time: {:.2f} seconds)'.format(0, timeit.default_timer() - start_time))
        df_dist = np.ones(df.shape) * -1
        df_dist = pd.DataFrame(data=df_dist, columns=col, index=idx)
        # print(len(idx_list), len(df_dist_list))

        for i in range(len(idx_list)):
            idx_local = df_dist_list[i].index
            col_local = df_dist_list[i].columns
            df_dist.loc[idx_local, col_local] = df_dist_list[i]

        filename_dis = filename.replace('.csv', '_distance.csv')

        print('CPU {:04d}: Now saving pandas DataFrame to {:s} (Time: {:.2f} seconds)'.format(0, filename_dis, timeit.default_timer() - start_time))
        df_dist = df_dist.astype(np.int8)
        df_dist.to_csv(filename_dis, sep=',', header=True, index=True)

        elapsed_time = timeit.default_timer() - start_time
        hour = math.floor(elapsed_time / 3600)
        minute = math.floor((elapsed_time - hour * 3600) / 60)
        second = elapsed_time - 3600 * hour - 60 * minute
        print('Time elapsed: {:d} hours {:d} minutes {:.2f} seconds'.format(hour, minute, second))

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
    Z = A.astype(np.float32).dot(A.astype(np.float32))
    B = np.array([[1 if i != j and ( A[i,j] == 1 or Z[i,j] > 0 ) else 0 for j in range(n)] for i in range(n)])
    T = APD_recursive(B)
    X = T.astype(np.float32).dot(A.astype(np.float32))
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
        Z = B.astype(np.float32).dot(B.astype(np.float32))
        B = np.array([[1 if i != j and ( B[i,j] == 1 or Z[i,j] > 0 ) else 0 for j in range(n)] for i in range(n)])
    T = B
    X = T.astype(np.float32).dot(A.astype(np.float32))
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
        # print(components[i], idx.shape, col.shape, type(idx), type(col))

    print('CPU {:04d}: Reconstructing adjacency pandas DataFrame (Time: {:.2f} seconds)'.format(0, timeit.default_timer() - start_time))
    df_orig = np.zeros(df.shape)
    df_orig = pd.DataFrame(data=df_orig, columns=col, index=idx)
    for i in range(len(idx_list)):
        df_orig.loc[idx_list[i], col_list[i]] = df.loc[idx_list[i], col_list[i]]
    filename_adj = filename.replace('.csv', '_adjacency.csv')
    print('CPU {:04d}: Now saving pandas DataFrame to {:s} (Time: {:.2f} seconds)'.format(0, filename_adj, timeit.default_timer() - start_time))
    df_orig = df_orig.astype(np.int8)
    df_orig.to_csv(filename_adj, sep=',', header=True, index=True)
    return (df, components, idx, col, idx_list, col_list)

def build_overlap_matrix_parallel(rows, filename='./overlap_matrix_parallel_block.csv'):
    start_time = timeit.default_timer()
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    status = MPI.Status()
    data = []
    n_indices = 0
    # print('Total CPUs: {:4d}, Rank {:4d}, rows: {:5d}, rows_per_cpu: {:5d}, row_start: {:5d}, cols_per_cpu: {:5d}, col_start: {:5d}'.format(size, rank, rows, rows_local, row_start, cols_local, col_start))

    if rank == 0:
        file_open = fu.FileOpen("data", "tweet-data.csv")
        tweet_spatial_analysis_config = cu.TweetSpatialAnalysisConfig("conf/tweet_spatial_analysis.ini")
        tdp = TweetDataPreProcessing(file_open, tweet_spatial_analysis_config)
        tdp.read_from_json("data/tweet_mean_all.json",
                           "data/tweets_median_working.json",
                           "data/tweets_median_non_working.json")
        df = tdp.tweet_data_working.df
        if rows == 0:
            rows, _ = df.shape
        mpi4py.get_config()
        np.__config__.show()
        tweet_data_working_overlap = np.zeros((rows, rows))
        indices = [(i, j) for i in range(rows) for j in range(rows) if j < i]
        n_indices = len(indices)
        n_indices = comm.bcast(n_indices, root=0)
        n_sent = 0
        for k in range(min(size - 1, n_indices)):
            i, j = indices[n_sent]
            data =  i, j, \
                    df['x'][i], df['y'][i], df['a'][i], df['b'][i], df['angle'][i], \
                    df['x'][j], df['y'][j], df['a'][j], df['b'][j], df['angle'][j]
            comm.send(data, dest=n_sent + 1, tag=n_sent)
            n_sent += 1

        for k in range(n_indices):
            overlap = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            src = status.Get_source()
            tag = status.Get_tag()

            i, j = indices[tag]
            tweet_data_working_overlap[i, j] = overlap
            if n_sent < n_indices:

                i, j = indices[n_sent]
                data =  i, j, \
                        df['x'][i], df['y'][i], df['a'][i], df['b'][i], df['angle'][i], \
                        df['x'][j], df['y'][j], df['a'][j], df['b'][j], df['angle'][j]
                comm.send(data, dest=src, tag=n_sent)
                n_sent += 1
            else:
                comm.send(indices[0], dest=src, tag=60000)
            # print('CPU {:04d}: Received data from CPU {:04d} (Time: {:.2f} seconds)'.format(rank,  src, timeit.default_timer() - start_time))

            # if ((k + 1) / n_indices * 100) % 10 < 0.01:
            if (k + 1) % 1000000 == 0:
                print('CPU {:04d}: {:6.2f}% accomplished (Time: {:.2f} seconds)'.format(rank,  (k + 1) / n_indices * 100, timeit.default_timer() - start_time))

        tweet_data_working_overlap = tweet_data_working_overlap + tweet_data_working_overlap.T
        print('CPU {:04d}: Merged received data to global numpy array (Time: {:.2f} seconds)'.format(rank, timeit.default_timer() - start_time))
        tweet_data_working_overlap = pd.DataFrame(data=tweet_data_working_overlap, columns=df['id'][0:rows], index=df['id'][0:rows])
        print('CPU {:04d}: Converted numpy array to pandas DataFrame (Time: {:.2f} seconds)'.format(rank, timeit.default_timer() - start_time))

        data = find_components(tweet_data_working_overlap, filename, start_time)
        df, components, idx, col, idx_list, col_list = data
        number_of_components = len(components)
        number_of_components = comm.bcast(number_of_components, root=0)
        print('CPU {:04d}: Found {:d} components (Time: {:.2f} seconds)'.format(rank, number_of_components, timeit.default_timer() - start_time))

        nsent = 0
        for i in range(min(size - 1, number_of_components)):
            comm.send(df.loc[idx_list[i],col_list[i]], dest=i + 1, tag=i)
            nsent += 1

        df_dist_list = []
        for i in range(number_of_components):
            D = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
            src = status.Get_source()
            tag = status.Get_tag()

            # print('CPU {:04d}: Received data from CPU {:04d} (Time: {:.2f} seconds)'.format(rank, src, timeit.default_timer() - start_time))
            df_dist_list.append(pd.DataFrame(data=D, columns=col_list[tag], index=idx_list[tag]))

            if nsent < number_of_components:
                comm.send(df.loc[idx_list[nsent],col_list[nsent]], dest=src, tag=nsent)
            else:
                comm.send([], dest=src, tag=number_of_components + 10)

            nsent += 1

        print('CPU {:04d}: Constructing distance pandas DataFrame (Time: {:.2f} seconds)'.format(0, timeit.default_timer() - start_time))
        df_dist = np.ones(df.shape) * -1
        df_dist = pd.DataFrame(data=df_dist, columns=col, index=idx)
        # print(len(idx_list), len(df_dist_list))

        for i in range(len(idx_list)):
            idx_local = df_dist_list[i].index
            col_local = df_dist_list[i].columns
            df_dist.loc[idx_local, col_local] = df_dist_list[i]

        filename_dis = filename.replace('.csv', '_distance.csv')

        print('CPU {:04d}: Now saving pandas DataFrame to {:s} (Time: {:.2f} seconds)'.format(0, filename_dis, timeit.default_timer() - start_time))
        df_dist = df_dist.astype(np.int8)
        df_dist.to_csv(filename_dis, sep=',', header=True, index=True)

        elapsed_time = timeit.default_timer() - start_time
        hour = math.floor(elapsed_time / 3600)
        minute = math.floor((elapsed_time - hour * 3600) / 60)
        second = elapsed_time - 3600 * hour - 60 * minute
        print('Time elapsed: {:d} hours {:d} minutes {:.2f} seconds'.format(hour, minute, second))

    if rank > 0:
        n_indices = comm.bcast(n_indices, root=0)
        if rank < n_indices:
            while True:
                data = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
                tag = status.Get_tag()
                if tag == 60000:
                    break
                i, j, x1, y1, a1, b1, t1, x2, y2, a2, b2, t2 = data
                if i == j:
                    overlap = 0
                elif j < i:
                    if a1 == 0 or b1 == 0 or a2 == 0 or b2 == 0:
                        overlap = -1
                    else:
                        overlap = au.are_two_ellipses_overlapping(x1, y1, a1, b1, t1, x2, y2, a2, b2, t2)
                comm.send(overlap, dest=0, tag=tag)
                # print('CPU {:04d}: Sent data to CPU {:04d}'.format(rank,  0))

        number_of_components = 0
        number_of_components = comm.bcast(number_of_components, root=0)
        if rank < number_of_components:
            while True:
                A = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
                tag = status.Get_tag()
                if tag == number_of_components + 10:
                    break
                m, n = A.shape
                if m > 1:
                    print('CPU {:04d}: Starting APD for {:d}-by-{:d} matrix (Time: {:.2f} seconds)'.format(rank, m, n, timeit.default_timer() - start_time))
                D = APD_recursive(np.array(A))
                comm.send(D, dest=0, tag=tag)


if __name__ == '__main__':
    dt = datetime.datetime.now()
    year = dt.year
    month = dt.month
    day= dt.day
    hour = dt.hour
    minute = dt.minute
    second = dt.second

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

    if (MPI.COMM_WORLD.Get_size() == 1):
        algorithm = 'serial'
    else:
        algorithm = 'parallel'

    name = 'overlap_matrix_{:s}_{:s}_{:s}_{:d}_{:02d}_{:02d}_{:02d}_{:02d}_{:02d}.csv'.format(node, partition, algorithm, year, month, day, hour, minute, second)

    if args.filename:
        filename = args.filename
    else:
        filename = '{:s}/{:s}'.format(loc, name)

    # rows, _ = tdp.tweet_data_working.df.shape
    if args.rows:
        rows = int(args.rows)
    else:
        rows = 0

    if (MPI.COMM_WORLD.Get_size() == 1):
        build_overlap_matrix(rows, filename)
    else:
        build_overlap_matrix_parallel(rows, filename)

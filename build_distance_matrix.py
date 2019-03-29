#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Libao Jin <jinlibao@outlook.com>
#
# Distributed under terms of the MIT license.

"""
Build distance matrix
"""

import numpy as np
import pandas as pd
import timeit


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

def assemble_index_columns(list_name):
    '''
    Reorder index/columns of pandas.DataFrame
    '''
    idx = list_name[0].index
    col = list_name[0].columns
    idx_list = [list_name[0].index]
    col_list = [list_name[0].columns]
    for i in range(1, len(list_name)):
        idx = idx.append(list_name[i].index)
        col = col.append(list_name[i].columns)
        idx_list.append(list_name[i].index)
        col_list.append(list_name[i].columns)

    return (idx, col, idx_list, col_list)

def build_distance_matrix(df, filename, start_time):
    # Data cleaning
    col = [idx for idx in df.iloc[0].index if df.iloc[0][idx] == -1]
    idx = [int(idx) for idx in col]
    df = df.drop(index=idx, columns=col)

    components = bfs_df(df)
    df_dist_list = []
    index_list = []
    columns_list = []

    for i in range(len(components)):
        idx = df.index[components[i]]
        col = df.columns[components[i]]
        D = APD(np.array(df.loc[idx,col]))
        df_dist_list.append(pd.DataFrame(data=D.astype(int), columns=col, index=idx))

    df_dist_list.sort(key=len, reverse=True)
    idx, col, idx_list, col_list = assemble_index_columns(df_dist_list)

    df_orig = np.zeros(df.shape)
    df_orig = pd.DataFrame(data=df_orig.astype(int), columns=col, index=idx)
    df_dist = np.ones(df.shape) * -1
    df_dist = pd.DataFrame(data=df_dist.astype(int), columns=col, index=idx)

    for i in range(len(idx_list)):
        df_dist.loc[idx_list[i], col_list[i]] = df_dist_list[i]
        df_orig.loc[idx_list[i], col_list[i]] = df.loc[idx_list[i], col_list[i]]

    filename_adj = filename.replace('.csv', '_adjacency.csv')
    filename_dis = filename.replace('.csv', '_distance.csv')

    print('CPU {:04d}: Now saving pandas DataFrame to {:s} (Time: {:.2f} seconds)'.format(0, filename_adj, timeit.default_timer() - start_time))
    df_orig.to_csv(filename_adj, sep=',', header=True, index=True)

    print('CPU {:04d}: Now saving pandas DataFrame to {:s} (Time: {:.2f} seconds)'.format(0, filename_dis, timeit.default_timer() - start_time))
    df_dist.to_csv(filename_dis, sep=',', header=True, index=True)


if __name__ == '__main__':
    build_distance_matrix()

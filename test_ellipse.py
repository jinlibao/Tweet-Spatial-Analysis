#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Libao Jin <jinlibao@outlook.com>
#
# Distributed under terms of the MIT license.

"""
Test ellipse
"""
from utils import analysis_utilities as au
from utils import configuration_utilities as cu
from utils import file_utilities as fu
from data_preprocessing import TweetDataPreProcessing
import build_overlap_matrix as bom
import pandas as pd
import numpy as np

def get_shortest_path_by_id(id_from, id_to, suc_file='~/Documents/data/twitter/csv/sample/tweets_median_working_successor_matrix.csv', id_file='~/Documents/data/twitter/csv/sample/tweets_median_working_adjacency_matrix_ordered_id.csv'):
    id_mat = pd.read_csv(id_file, header=None)
    S = pd.read_csv(suc_file, header=None)

    idx_from = id_mat[id_mat.iloc[:, 1] == id_from].index.values[0]
    idx_to = id_mat[id_mat.iloc[:, 1] == id_to].index.values[0]

    path_nodes = [id_mat.iloc[idx_from, 1]]
    while idx_from != idx_to:
        idx_from = S.iloc[idx_from, idx_to]
        path_nodes.append(id_mat.iloc[idx_from, 1])
    print(path_nodes)
    return path_nodes

def get_shortest_path_by_idx(idx_from, idx_to, suc_file='~/Documents/data/twitter/csv/sample/tweets_median_working_successor_matrix.csv', id_file='~/Documents/data/twitter/csv/sample/tweets_median_working_adjacency_matrix_ordered_id.csv'):
    id_mat = pd.read_csv(id_file, header=None)
    S = pd.read_csv(suc_file, header=None)

    path_nodes = [id_mat.iloc[idx_from, 1]]
    while idx_from != idx_to:
        idx_from = S.iloc[idx_from, idx_to]
        path_nodes.append(id_mat.iloc[idx_from, 1])
    print(path_nodes)
    return path_nodes

def check_shortest_path_by_id(ids, raw_file='~/Documents/data/twitter/csv/sample/tweets_median_working.csv'):
    df = pd.read_csv(raw_file, header=None, names=['x', 'y', 'a', 'b', 'angle', 'id', 'area'])
    for i in ids:
        id_from, id_to = i
        path_nodes = get_shortest_path_by_id(id_from, id_to)
        data = [list(df[df['id'] == node][['x', 'y', 'a', 'b', 'angle', 'id']].values[0]) for node in path_nodes]
        au.plot_shortest_path(data, '{:d}_{:d}.pdf'.format(path_nodes[0], path_nodes[-1]))

def check_shortest_path_by_idx(idx, raw_file='~/Documents/data/twitter/csv/sample/tweets_median_working.csv'):
    df = pd.read_csv(raw_file, header=None, names=['x', 'y', 'a', 'b', 'angle', 'id', 'area'])
    for i in idx:
        idx_from, idx_to = i
        path_nodes = get_shortest_path_by_idx(idx_from, idx_to)
        data = [list(df[df['id'] == node][['x', 'y', 'a', 'b', 'angle', 'id']].values[0]) for node in path_nodes]
        # au.plot_shortest_path(data, './plot/{:d}_{:d}.pdf'.format(path_nodes[0], path_nodes[-1]))
        au.plot_shortest_path(data, './plot/{:d}_{:d}.pdf'.format(idx_from, idx_to))

def check_ellipse_by_id(id_pair, raw_file='~/Documents/data/twitter/csv/sample/tweets_median_working.csv'):
    df = pd.read_csv(raw_file, header=None, names=['x', 'y', 'a', 'b', 'angle', 'id', 'area'])
    data = [list(df[df['id'] == i][['x', 'y', 'a', 'b', 'angle', 'id']].values[0]) for i in id_pair]
    print(data)
    x0, y0, a0, b0, phi0, id0 = data[0]
    x1, y1, a1, b1, phi1, id1 = data[1]
    au.plot_ellipses(x0, y0, a0, b0, phi0, x1, y1, a1, b1, phi1, '{:d}_{:d}.pdf'.format(id_pair[0], id_pair[-1]))
    print(au.are_two_ellipses_overlapping(x0, y0, a0, b0, phi0, x1, y1, a1, b1, phi1))

def test_ellipse():
    A = [
        # [0, 0, 2, 1, 0, 0, 0, 3, 2, 0],
        # [2, 1, 2, 1, 0, -1, 2, 3, 2, 0.5],
        # [2, 1, 2, 1, 0.5, -1, 2, 3, 1.5, 0.5],
        # [2, 1, 2, 1, 0.5, -1, 2, 3, 1.5, 1.6],
        # [-8130485.626946354,4991209.679861832,117783.0664383206,3364.3910579424,-1.523861468, -8199215.743300232,5012710.176886535,77418.8445014143,29759.4265258827,0.43077731700000005],
        # [-8623270.760212664,4566528.186014481,31555.9614987463,3384.7972760126,-0.093125154, -8633386.872184316,4500582.478387021,71903.7423332023,32923.3000556941,1.328446492]
        [-8571408.115967816,4712938.215549811,4662.7670058175,2478.4635687348,1.5233426159999999, -8569898.603635151,4707530.330853039,4427.4039101461,3468.2101122412,0.699779302]
    ]
    for i in range(len(A)):
        x1, y1, a1, b1, phi1, x2, y2, a2, b2, phi2 = A[i]
        print(au.are_two_ellipses_overlapping(x1, y1, a1, b1, phi1, x2, y2, a2, b2, phi2))
        au.plot_ellipses(x1, y1, a1, b1, phi1, x2, y2, a2, b2, phi2)

def generate_csv(raw_files, csv_files):
    all_file, working_file, non_working_file = raw_files
    file_open = fu.FileOpen("data", "tweet-data.csv")
    tweet_spatial_analysis_config = cu.TweetSpatialAnalysisConfig("conf/tweet_spatial_analysis.ini")
    tdp = TweetDataPreProcessing(file_open, tweet_spatial_analysis_config)
    tdp.read_from_json(all_file, working_file, non_working_file)
    df = []
    df.append(tdp.tweet_data_working.df[['x', 'y', 'a', 'b', 'angle', 'id', 'area']])
    df.append(tdp.tweet_data_non_working.df[['x', 'y', 'a', 'b', 'angle', 'id', 'area']])
    df.append(tdp.tweet_data_all.df[['x', 'y', 'a', 'b', 'angle', 'id', 'area']])
    for i in range(len(df)):
        df[i].to_csv(csv_files[i], sep=',', header=False, index=False)

def test_generate_csv():
    src = 'sample'
    des = 'sample'
    # src = 'data'
    # des = 'data'
    names  = ["tweets_mean_all.json", "tweets_median_working.json", "tweets_median_non_working.json"]
    raw_files = [src + '/' + name for name in names]
    csv_files = [des + '/' + name.replace('json', 'csv') for name in names]
    generate_csv(raw_files, csv_files)

def test_check_ellipse():
    ids = [[2355845772, 1587067812], [2355845772, 14289663], [2573296611, 15067366], [25620815, 700476792], [419905564, 700476792], [15067366, 700476792], [190897931, 25620815],
            [259795124, 17179749], [33252934, 2573296611], [741621378, 15020935], [30380284, 15020935]]
    check_shortest_path_by_id(ids)

    ids = [[2573296611, 15067366]]
    for id_pair in ids:
        check_ellipse_by_id(id_pair)

    idx = [[i, j] for j in range(875,918) for i in range(875,j)]
    check_shortest_path_by_idx(idx)


if __name__ == '__main__':
    test_ellipse()
    # test_generate_csv()
    # test_check_ellipse()



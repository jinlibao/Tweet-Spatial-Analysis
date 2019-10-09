#! /usr/bin/env python3
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2019 Libao Jin <jinlibao@outlook.com>
#
# Distributed under terms of the MIT license.
"""
Generate random matrix
"""

import numpy as np
import pandas as pd


def create_matrix(n):
    a = np.random.randint(0, 10, (n, n))
    a = a + a.T
    b = np.random.randint(0, 10, (n, n))
    b = b + b.T
    c = a.dot(b)
    aa = a.dot(a)
    A = pd.DataFrame(data=a)
    B = pd.DataFrame(data=b)
    C = pd.DataFrame(data=c)
    AA = pd.DataFrame(data=aa)
    A.to_csv("data_{:d}_A.csv".format(n), index=False, header=False)
    B.to_csv("data_{:d}_B.csv".format(n), index=False, header=False)
    C.to_csv("data_{:d}_C.csv".format(n), index=False, header=False)
    AA.to_csv("data_{:d}_AA.csv".format(n), index=False, header=False)
    print(A)
    print(B)
    print(C)


if __name__ == '__main__':
    create_matrix(10)
    create_matrix(100)
    create_matrix(1000)
    create_matrix(2000)

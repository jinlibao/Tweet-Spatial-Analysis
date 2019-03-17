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

def test_ellipse():
    # x1, y1, a1, b1, phi1, x2, y2, a2, b2, phi2 = 0, 0, 2, 1, 0, 0, 0, 3, 2, 0
    # x1, y1, a1, b1, phi1, x2, y2, a2, b2, phi2 = 2, 1, 2, 1, 0, -1, 2, 3, 2, 0.5
    # x1, y1, a1, b1, phi1, x2, y2, a2, b2, phi2 = 2, 1, 2, 1, 0.5, -1, 2, 3, 1.5, 0.5
    x1, y1, a1, b1, phi1, x2, y2, a2, b2, phi2 = 2, 1, 2, 1, 0.5, -1, 2, 3, 1.5, 1.6
    print(au.are_two_ellipses_overlapping(x1, y1, a1, b1, phi1, x2, y2, a2, b2, phi2))
    au.plot_ellipses(x1, y1, a1, b1, phi1, x2, y2, a2, b2, phi2)


if __name__ == '__main__':
    test_ellipse()

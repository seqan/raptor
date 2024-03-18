#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

import numpy as np
import pandas as pd

import os
import re

PATTERN_SIZE=250
FREQUENCIES_DIR = '<path>/{}_p2'.format(PATTERN_SIZE)

def get_wk_list(path = FREQUENCIES_DIR):
    result = []
    for entry in os.listdir(path):
        match = re.match('(\d+)_(\d+)_(\d+\w).out.minimiser_counts$', entry)
        if match is not None:
            result.append( (match.group(1), match.group(2), match.group(3)) )
    return result

def generate_table():
    data = []
    params = get_wk_list()
    for (window_size, kmer_size, ibf_size) in params:
        column_name = 'w{} k{}'.format(window_size, kmer_size)
        path_to_threshold = os.path.join(FREQUENCIES_DIR, '{}_{}_{}.out.minimiser_counts'.format(window_size, kmer_size, ibf_size))
        data.append(pd.read_csv(path_to_threshold, delimiter='\t', header=0, index_col=False, usecols=[1], names=[column_name]))
    output_file = os.path.join(FREQUENCIES_DIR, 'p{}_w{}_k{}.frequencies'.format(PATTERN_SIZE, window_size, kmer_size))
    pd.concat(data, axis=1, copy=False).to_csv(output_file, float_format='%d')

generate_table()

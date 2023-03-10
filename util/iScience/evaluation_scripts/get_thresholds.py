#!/usr/bin/env python3

# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

import numpy as np
import pandas as pd

import os
import re

THRESHOLD_DIR = '<path>/250_e2'

def get_wk_list(path = THRESHOLD_DIR):
    result = []
    for entry in os.listdir(path):
        match = re.match('text_p(\d+)_w(\d+)_k(\d+)_e(\d+)_tau0.txt$', entry)
        if match is not None:
            result.append( (match.group(1), match.group(2), match.group(3), match.group(4)) )
    return result

def generate_table():
    data = []
    params = get_wk_list()
    for (pattern_size, window_size, kmer_size, error_count) in params:
        column_name = 'w{} k{}'.format(window_size, kmer_size)
        path_to_threshold = os.path.join(THRESHOLD_DIR, 'text_p{}_w{}_k{}_e{}_tau0.txt'.format(pattern_size, window_size, kmer_size, error_count))
        data.append(pd.read_csv(path_to_threshold, delimiter='\t', header=None, index_col=False, usecols=[1], names=[column_name]))
    (pattern_size, _, _, error_count) = params[0]
    output_file = os.path.join(THRESHOLD_DIR, 'p{}_e{}_tau09999.threshold'.format(pattern_size, error_count))
    pd.concat(data, axis=1, copy=False).to_csv(output_file, float_format='%d')

generate_table()

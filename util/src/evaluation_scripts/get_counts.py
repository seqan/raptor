#!/usr/bin/env python3

# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

import numpy as np
import pandas as pd
from math import ceil

import os
import re

COUNTS_DIR = '<output_of_counting>' # directory containing w_k.counts, e.g. 23_19.counts

def get_wk_list(path = COUNTS_DIR):
    result = []
    for entry in os.listdir(path):
        match = re.match('(\d+)_(\d+).counts$', entry)
        if match is not None:
            result.append( (match.group(1), match.group(2)) )
    return sorted(result)

def process_file(path):
    with open(path, 'r') as file:
        total_text_size = int(file.readline().strip())
        bin_sizes = []
        for line in file:
            bin_sizes.append(int(line.strip()))
    return (total_text_size, max(bin_sizes), ceil(sum(bin_sizes)/len(bin_sizes)), sum(bin_sizes))

def generate_table():
    data = []
    params = get_wk_list()
    for (window_size, kmer_size) in params:
        (total_text, max, avg, sum) = process_file(os.path.join(COUNTS_DIR, '{}_{}.counts'.format(window_size, kmer_size)))
        data.append([total_text, max, avg, sum])
    df = pd.DataFrame(data, columns = ['text size', 'max bin', 'avg bin', 'sum bin'])
    df.index = ['w{} k{}'.format(window_size, kmer_size) for [window_size, kmer_size] in params]
    output_file = os.path.join(COUNTS_DIR, 'counts.csv')
    df.to_csv(output_file)

generate_table()

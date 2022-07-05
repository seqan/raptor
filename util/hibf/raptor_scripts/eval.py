#!/usr/bin/env python3

# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

import numpy as np
import pandas as pd

import os
import re
import time

BINS=65536
READ_COUNT=2**20
REPEATS=10
INPUT_BASE = '/example/runs'

def process_output(path):
    with open(path) as f:
        tp = 0
        fp = 0
        fn = 0
        bin_to_ub = {}
        for line in f:
            if line[0] == '#':
                try:
                    [ub, bin] = line.strip().split('\t')
                    ub = int(ub[1:])
                    bin = int(re.match('bin_(\d+)\..*$', os.path.basename(bin))[1])
                    bin_to_ub[bin] = ub
                except ValueError:
                    pass
                continue
            try:
                [x, y] = line.strip().split('\t')
            except ValueError:
                fn += 1
                continue
            [read_id, bins] = [int(x), [int(e) for e in y.split(',') if e != '']]
            true_id = bin_to_ub[(read_id % READ_COUNT) // (READ_COUNT // BINS)]
            if true_id in bins:
                tp += 1
                if len(bins) != 1:
                    fp += len(bins) - 1
            else:
                fn += 1
                fp += len(bins)
        return [int(tp/REPEATS),int(fp/REPEATS),int(fn/REPEATS)]

def process_time(path):
    time = []
    mem = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i == 4:
                time = line.strip().split(' ')[-1].split(':')
                if len(time) == 3:
                    time = '{}:{:02d}'.format(int(time[0]) * 60 + int(time[1]), int(time[2]))
                elif len(time) == 2:
                    if int(time[0]) == 0:
                        time = '{}:{}'.format(time[0], time[1])
                    else:
                        time = '{}:{:02d}'.format(time[0], round(float(time[1])))
            if i == 9:
                mem = round(int(line.strip().split(' ')[-1]) / 1024)
    return [time, mem]

def process_internal_time(path):
    with open(path) as f:
        line = f.readline()
        line = f.readline()
        [ibf_time, read_time, compute_time] = line.strip().split('\t')
        return [ibf_time, read_time, compute_time]

def get_param_list(path):
    result = []
    for entry in os.listdir(path):
        match = re.match('(\d+)_(\d+)_(\w+).out$', entry)
        if match is not None:
            result.append( (match.group(1), match.group(2), match.group(3)) )
    result.sort()
    return result

def generate_table(path):
    data = []
    params = get_param_list(path)
    format_string = 'Processing file {{:>{}}} of {}...'.format(len(str(len(params))), len(params))
    for i, (window_size, kmer_size, tmax) in enumerate(params):
        print(format_string.format(i + 1), end='', flush=True)
        path_to_build_time = os.path.join(path, '{}_{}_{}_build.time'.format(window_size, kmer_size, tmax))
        path_to_query_time = os.path.join(path, '{}_{}_{}_query.time'.format(window_size, kmer_size, tmax))
        path_to_output = os.path.join(path, '{}_{}_{}.out'.format(window_size, kmer_size, tmax))
        path_to_internal_time = os.path.join(path, '{}_{}_{}.out.time'.format(window_size, kmer_size, tmax))

        # count_time = process_time(path_to_count_time)
        # layout_time =
        [build_time, build_ram] = process_time(path_to_build_time)
        [query_time, query_ram] = process_time(path_to_query_time)
        [query_ibf_time, query_reads_time, query_compute_time] = process_internal_time(path_to_internal_time)
        [tp, fp, fn] = [-1,-1,-1]#process_output(path_to_output)

        data.append([build_time,
                     build_ram,
                     query_time,
                     query_ibf_time,
                     query_reads_time,
                     query_compute_time,
                     query_ram,
                     fp,
                     fn])
        print('Done', flush=True)
    df = pd.DataFrame(data, columns = [('Construct', 'Time [MM:SS]'),
                                       ('Construct', 'RAM [MiB]'),
                                       ('Search', 'Overall [MM:SS.ss]'),
                                       ('Search', 'IBF I/O [SS.ss]'),
                                       ('Search', 'Reads I/O [SS.ss]'),
                                       ('Search', 'Compute [SS.ss]'),
                                       ('Search', 'RAM [MiB]'),
                                       ('Search', 'FP'),
                                       ('Search', 'FN')])
    df.index = [str(x).replace(',', ';').replace("'", '') for x in params]
    df.columns = pd.MultiIndex.from_tuples(df.columns, names=['','(w; k; size)'])
    return df

# input_path = os.path.join(INPUT_BASE, str(BINS))
# if not os.path.exists(input_path):
#     raise OSError("{} does not exist.".format(input_path))
# output_path = os.path.join(INPUT_BASE, str(BINS))
# os.makedirs(output_path, exist_ok=True)

# path_to_csv = os.path.join(output_path, 'table.csv')
# df = generate_table('/example/runs/raptor/65536/0.05')
df = generate_table('/example/runs/raptor/raptor_bench/65536/0.015')
print(df)
# df.to_csv(path_to_csv)

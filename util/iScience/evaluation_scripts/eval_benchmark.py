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
import time

BINS=1024
READ_COUNT=1048576
INPUT_BASE = '/home/infri/develop/raptor/build/util/test/results' # E.g., /dev/shm/username; should contain BINS subdirectory.
OUTPUT_BASE = '/home/infri/develop/raptor/build/util/test/results' # Will create BINS subdirectory.
EVAL_ENERGY=False # Energy consumption was measured in the benchmark.

input_path = os.path.join(INPUT_BASE, str(BINS))
if not os.path.exists(input_path):
    raise OSError("{} does not exist.".format(input_path))
output_path = os.path.join(INPUT_BASE, str(BINS))
os.makedirs(output_path, exist_ok=True)

def process_output(path):
    with open(path) as f:
        tp = 0
        fp = 0
        fn = 0
        for line in f:
            if line[0] == '#':
                continue
            try:
                [x, y] = line.strip().split('\t')
            except ValueError:
                fn += 1
                continue
            [read_id, bins] = [int(x), [int(e) for e in y.split(',') if e != '']]
            true_id = (read_id % READ_COUNT) // (READ_COUNT // BINS)
            if true_id in bins:
                tp += 1
                if len(bins) != 1:
                    fp += len(bins) - 1
            else:
                fn += 1
                fp += len(bins)
        return [tp,fp,fn]

def process_time_log(path):
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

def process_energy(path):
    if not EVAL_ENERGY:
        return [-1, -1]
    with open(path) as f:
        line = f.readline()
        line = f.readline()
        line = f.readline()
        line = f.readline()
        line = f.readline()
        line = f.readline()
        energy_pkg = line.strip().split(' ')[0]
        line = f.readline()
        energy_ram = line.strip().split(' ')[0]
        return [energy_pkg, energy_ram]

def get_param_list(path = output_path):
    result = []
    for entry in os.listdir(path):
        match = re.match('(\d+)_(\d+)_(\d+\w).out$', entry)
        if match is not None:
            result.append( (match.group(1), match.group(2), match.group(3)) )
    result.sort(key=lambda entry: int(entry[2][:-1]))
    return result

def generate_table_unpartitioned():
    data = []
    params = get_param_list()
    format_string = 'Processing file {{:>{}}} of {}...'.format(len(str(len(params))), len(params))
    for i, (window_size, kmer_size, ibf_size) in enumerate(params):
        print(format_string.format(i + 1), end='', flush=True)
        path_to_build_log = os.path.join(output_path, '{}_{}_{}_build.log'.format(window_size, kmer_size, ibf_size))
        path_to_build_perf = os.path.join(input_path, '{}_{}_{}_build.perf'.format(window_size, kmer_size, ibf_size))
        path_to_query_log = os.path.join(output_path, '{}_{}_{}_query.log'.format(window_size, kmer_size, ibf_size))
        path_to_query_perf = os.path.join(input_path, '{}_{}_{}_query.perf'.format(window_size, kmer_size, ibf_size))
        path_to_output = os.path.join(output_path, '{}_{}_{}.out'.format(window_size, kmer_size, ibf_size))
        path_to_internal_time = os.path.join(output_path, '{}_{}_{}.out.time'.format(window_size, kmer_size, ibf_size))

        [build_time, build_ram] = process_time_log(path_to_build_log)
        [build_energy_pkg, build_energy_ram] = process_energy(path_to_build_perf)
        [query_time, query_ram] = process_time_log(path_to_query_log)
        [query_energy_pkg, query_energy_ram] = process_energy(path_to_query_perf)
        [query_ibf_time, query_reads_time, query_compute_time] = process_internal_time(path_to_internal_time)
        [tp, fp, fn] = process_output(path_to_output)

        data.append([build_time,
                     build_ram,
                     build_energy_pkg,
                     build_energy_ram,
                     query_time,
                     query_ibf_time,
                     query_reads_time,
                     query_compute_time,
                     query_ram,
                     query_energy_pkg,
                     query_energy_ram,
                     fp,
                     fn])
        print('Done', flush=True)
    df = pd.DataFrame(data, columns = [('Construct', 'Time [MM:SS]'),
                                       ('Construct', 'RAM [MiB]'),
                                       ('Construct', 'energy-pkg [J]'),
                                       ('Construct', 'energy-ram [J]'),
                                       ('Search', 'Overall [MM:SS.ss]'),
                                       ('Search', 'IBF I/O [SS.ss]'),
                                       ('Search', 'Reads I/O [SS.ss]'),
                                       ('Search', 'Compute [SS.ss]'),
                                       ('Search', 'RAM [MiB]'),
                                       ('Search', 'energy-pkg [J]'),
                                       ('Search', 'energy-ram [J]'),
                                       ('Search', 'FP'),
                                       ('Search', 'FN')])
    df.index = [str(x).replace(',', ';').replace("'", '') for x in params]
    df.columns = pd.MultiIndex.from_tuples(df.columns, names=['','(w; k; size)'])
    if not EVAL_ENERGY:
        df = df.drop([('Construct', 'energy-pkg [J]'), ('Construct', 'energy-ram [J]'),
                      ('Search', 'energy-pkg [J]'), ('Search', 'energy-ram [J]')], axis=1)
    return df

path_to_csv = os.path.join(output_path, 'table.csv')
df = generate_table_unpartitioned()
df.to_csv(path_to_csv)

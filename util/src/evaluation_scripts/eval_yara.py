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

from operator import itemgetter, attrgetter

# Adjust parameters in generate_table (bin and read count)
# stores results in `OUTPATH/table.csv`
# This directory contains all '.log' and '.time' files for a run
OUT_PATH = '<output directory with bin_number, e.g. /dev/shm/username/1024'
# Index is only build once, provide path to time log
path_to_build_fm_time = os.path.join(OUT_PATH, '19_32G_build_fm.time')

# Original DREAM-Yara does not use minimiser -> False
# New DREAM-Yara uses minimiser -> True
MINIMISERS=False

def process_time_log(path):
    time = []
    mem = []
    with open(path) as f:
        for i, line in enumerate(f):
            if i == 4:
                time = line.strip().split(' ')[-1].split(':')
                if len(time) == 3:
                    time = '{}:{}:{:02d}'.format(int(time[0]), int(time[1]), int(time[2]))
                elif len(time) == 2:
                    if int(time[0]) == 0:
                        time = '{}:{}'.format(time[0], time[1])
                    else:
                        time = '{}:{:02d}'.format(time[0], round(float(time[1])))
            if i == 9:
                mem = round(int(line.strip().split(' ')[-1]) / 1024)
    return [time, mem]

def process_mapper_log(path):
    with open(path) as f:
        for line in f:
            if line.startswith('Filter loading time'):
                components = line.split('\t')
                ibf_load_time = round(float(components[2].rstrip(' sec')), 1)
            if line.startswith('Reads filtering time'):
                components = line.split('\t')
                ibf_filter_time = round(float(components[2].rstrip(' sec')), 1)
            if line.startswith('Total reads'):
                components = line.split('\t')
                total_reads = int(float(components[3].strip()))
            if line.startswith('Mapped reads'):
                components = line.split('\t')
                mapped_reads = int(float(components[3].strip()))
            if line.startswith('Avg reads per bin'):
                components = line.split('\t')
                avg_per_bin = int(float(components[2].strip()))
    return [ibf_load_time, ibf_filter_time, total_reads - mapped_reads, avg_per_bin]

def get_param_list(path = OUT_PATH):
    result = []
    for entry in os.listdir(path):
        if MINIMISERS:
            match = re.match('(\d+)_(\d+)_(\d+)G_mapper_(\d+).time$', entry)
        else:
            match = re.match('(\d+)_(\d+)G_mapper_(\d+).time$', entry)
        if match is not None:
            if MINIMISERS:
                result.append( (int(match.group(1)), int(match.group(2)), int(match.group(3)), int(match.group(4))) )
            else:
                result.append( (int(match.group(1)), int(match.group(2)), int(match.group(3))) )
    return result

# (how many bins: int), (how many reads: int)
def generate_table(bin_count=1024, read_count=1048576):
    data = []
    if MINIMISERS:
        params = sorted(get_param_list(), key=itemgetter(0,3,2))
    else:
        params = sorted(get_param_list(), key=itemgetter(0,2,1))
    format_string = 'Processing file {{:>{}}} of {}...'.format(len(str(len(params))), len(params))

    [indexer_time, indexer_ram] = process_time_log(path_to_build_fm_time)

    for i, tuple in enumerate(params):
        print(format_string.format(i + 1), end='', flush=True)
        # path_to_build_ibf_log = os.path.join(OUT_PATH, '{}_{}_{}G_build_ibf.log'.format(window_size, kmer_size, ibf_size))
        if MINIMISERS:
            window_size, kmer_size, ibf_size, read_size = tuple
            path_to_build_ibf_time = os.path.join(OUT_PATH, '{}_{}_{}G_build_ibf.time'.format(window_size, kmer_size, ibf_size))
            path_to_mapper_time = os.path.join(OUT_PATH, '{}_{}_{}G_mapper_{}.time'.format(window_size, kmer_size, ibf_size, read_size))
            path_to_mapper_log = os.path.join(OUT_PATH, '{}_{}_{}G_mapper_{}.log'.format(window_size, kmer_size, ibf_size, read_size))
        else:
            kmer_size, ibf_size, read_size = tuple
            path_to_build_ibf_time = os.path.join(OUT_PATH, '{}_{}G_build_ibf.time'.format(kmer_size, ibf_size))
            path_to_mapper_time = os.path.join(OUT_PATH, '{}_{}G_mapper_{}.time'.format(kmer_size, ibf_size, read_size))
            path_to_mapper_log = os.path.join(OUT_PATH, '{}_{}G_mapper_{}.log'.format(kmer_size, ibf_size, read_size))

        [build_ibf_time, build_ibf_ram] = process_time_log(path_to_build_ibf_time)
        [mapper_time, mapper_ram] = process_time_log(path_to_mapper_time)

        [mapper_ibf_load_time, mapper_ibf_filter_time, mapper_unmapped_reads, mapper_avg_per_bin] = process_mapper_log(path_to_mapper_log)

        data.append([build_ibf_time, build_ibf_ram,
                     indexer_time, indexer_ram,
                     mapper_time, mapper_ibf_load_time, mapper_ibf_filter_time, mapper_ram, mapper_unmapped_reads, mapper_avg_per_bin])
        print('Done', flush=True)
    df = pd.DataFrame(data, columns = [('IBF', 'Time'),  ('IBF', 'RAM'),
                                       ('Indexer','Time'),  ('Indexer','RAM'),
                                       ('Mapper','Time'),  ('Mapper','IBF Load Time'), ('Mapper','IBF Filter Time'),  ('Mapper','RAM'), ('Mapper','Unmapped Reads'), ('Mapper', 'Processed Reads per Bin')])
    df.index = [str(x).replace(',', ';').replace("'", '').replace('m', '') for x in params]
    if MINIMISERS:
        df.columns = pd.MultiIndex.from_tuples(df.columns, names=['','w; k; IBF size; read length'])
    else:
        df.columns = pd.MultiIndex.from_tuples(df.columns, names=['','k; IBF size; read length'])
    return df

path_to_csv = os.path.join(OUT_PATH, 'table.csv')
df = generate_table()
df.to_csv(path_to_csv)

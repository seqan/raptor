#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

#
# Usage ram_usage.py <input_file> <output_file>
#
# Computes a table with RAM-Usage from a file containing output of `time -v`.
import argparse
import os
import pandas

parser = argparse.ArgumentParser(description='Parse time and memory consumption of compiling.', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('input', type=str, help='File containing all outputs of `time -v`.')
parser.add_argument('output', type=str, help='File to write output to. (CSV)')
arguments = parser.parse_args()

file_names = []
ram_usages = []
run_times = []

with open(arguments.input, 'r') as input_file:
    parsing_ram_usage = False
    for line_number, line in enumerate(input_file):
        if line_number % 23 == 0:
            index_of_unit = line.rfind('-c')
            if index_of_unit != - 1:
                parsing_ram_usage = True
                file_names.append(line[index_of_unit:][:-2].split('/')[-1])
            else:
                parsing_ram_usage = False
        if parsing_ram_usage and ((line_number - 9) % 23) == 0:
            ram_usages.append(int(line.split(' ')[-1]) // 1024)
        if parsing_ram_usage and ((line_number - 4) % 23) == 0:
            run_times.append(line.strip().split(' ')[-1].lstrip('0:'))

with open(arguments.output, 'w') as output_file:
    df = pandas.DataFrame({'File' : file_names, 'RAM in MiB' : ram_usages, 'Time in s' : run_times})
    df = df.sort_values(by=['File'], ascending=True)
    df.to_csv(output_file, index=False)

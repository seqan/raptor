#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -e

BINARY_DIR="<path to built binaries>"
OUT_DIR="<output path>" # will create one file for each w/k combination
BIN_DIR="<bin path>" # output directory of simulation. the directory that contains the BIN_NUMBER directory
BIN_NUMBER=1024
THREADS=4 # note that multiple (5) processes are spawned (i.e., 5*4 threads), see line 27

output_dir=$OUT_DIR/$BIN_NUMBER
bin_dir=$BIN_DIR/$BIN_NUMBER/bins

mkdir -p $output_dir

do_task () {
    echo "Counting $w $k"
    $BINARY_DIR/count_minimiser \
        --window $w \
        --kmer $k \
        --threads $THREADS \
        --output $output_dir/$w\_$k.counts \
        $(seq -f "$bin_dir/bin_%0${#BIN_NUMBER}.0f.fasta" 0 1 $((BIN_NUMBER-1)))
}

pidlist=""

for w in $(seq 23 2 32 && seq 32 2 80) # Adjust for different w/k
do
    for k in 16 17 18 19 20
    do
        do_task & pidlist="$pidlist $!"
    done
    for job in $pidlist
    do
        wait $job
    done
done

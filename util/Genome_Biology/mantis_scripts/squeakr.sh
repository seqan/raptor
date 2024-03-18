#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -e
SCRIPT_ROOT=$(dirname $(readlink -f $0))
source $SCRIPT_ROOT/variables.sh

run_squeakr () {
    $SQUEAKR_BINARY count \
        -e \
        -k $K \
        -c 1 \
        -n \
        -t 1 \
        -o $SQUEAKR_DIRECTORY/$(basename ${FASTQ} .fastq).squeakr \
        ${FASTQ}
}

for start in $(seq 0 32 $((BIN_NUMBER-1))); do
    pidlist=""
    for FASTQ in $(seq -f "$WORKING_DIRECTORY/bins/bin_%0${#BIN_NUMBER}.0f.fastq" $start 1 $((start+31))); do
        run_squeakr & pidlist="$pidlist $!"
    done
    for job in $pidlist; do
        wait $job
    done
done


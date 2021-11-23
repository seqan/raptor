#!/usr/bin/env bash
# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

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
    for FASTQ in $(seq -f "$WORKING_DIRECTORY/bins/bin_%0${#BIN_NUMBER}g.fastq" $start 1 $((start+31))); do
        run_squeakr & pidlist="$pidlist $!"
    done
    for job in $pidlist; do
        wait $job
    done
done


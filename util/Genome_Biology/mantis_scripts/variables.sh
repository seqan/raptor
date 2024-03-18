#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -e

READ_LENGTH=100
ERRORS=2
K=20
THREADS=32
BIN_NUMBER=65536
SQUEAKR_BINARY="/example/build/squeakr/squeakr"
MANTIS_BINARY="/example/build/mantis/bin/mantis"
HELPER_BINARY="/example/build/helper/bin/fasta_to_fastq"
INPUT_DIR="/example/big_dataset"
BENCHMARK_DIR="/example/mantis_bench"

WORKING_DIRECTORY=$BENCHMARK_DIR/$BIN_NUMBER
mkdir -p $WORKING_DIRECTORY/bins/
mkdir -p $WORKING_DIRECTORY/reads/

SQUEAKR_DIRECTORY=$WORKING_DIRECTORY/squeakr
mkdir -p $SQUEAKR_DIRECTORY

MANTIS_INDEX=$WORKING_DIRECTORY/mantis/
mkdir -p $MANTIS_INDEX

READ_FILE=$WORKING_DIRECTORY/reads/all_10.fastq
QUERY_OUT=$WORKING_DIRECTORY/result.txt

copy_time=$WORKING_DIRECTORY/R$READ_LENGTH\_K$K\_copy.time
copy_log=$WORKING_DIRECTORY/R$READ_LENGTH\_K$K\_copy.log
squeakr_time=$WORKING_DIRECTORY/R$READ_LENGTH\_K$K\_squeakr.time
squeakr_log=$WORKING_DIRECTORY/R$READ_LENGTH\_K$K\_squeakr.log
mantis_build_time=$WORKING_DIRECTORY/R$READ_LENGTH\_K$K\_mantis_build.time
mantis_build_log=$WORKING_DIRECTORY/R$READ_LENGTH\_K$K\_mantis_build.log
mantis_mst_time=$WORKING_DIRECTORY/R$READ_LENGTH\_K$K\_mantis_mst.time
mantis_mst_log=$WORKING_DIRECTORY/R$READ_LENGTH\_K$K\_mantis_mst.log
mantis_query_time=$WORKING_DIRECTORY/R$READ_LENGTH\_K$K\_mantis_query.time
mantis_query_log=$WORKING_DIRECTORY/R$READ_LENGTH\_K$K\_mantis_query.log

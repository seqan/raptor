#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -Eeuo pipefail

echo "This script is no longer supported. Please use the new simulation scripts in the util/simulation directory."
exit 1

BINARY_DIR="<path to built binaries>" # Dir containing "mason_genome", "split_sequence", etc.
OUT_DIR="<output path>" # Where simulated data should be stored
LENGTH=536870912 # 512MiB
SEED=42 # was 20181406 before, but was hardcoded to 42 in seqan

# Differences:
# * No haplotype
# * Only reads of length 250
# * No all_10.fastq

# Requirements
# split_sequence.cpp:105
# std::string filename = "bin_" + padded_parts + ".fa"
# to
# std::string filename = "bin_" + padded_parts + ".fasta"

for EXPONENT in $(seq 10 1 20); do
    BIN_NUMBER=$(bc <<< "2^${EXPONENT}")
    ERRORS=2
    READ_LENGTH=250
    READ_COUNT=1048576
    HAPLOTYPE_COUNT=1

    output_dir=$OUT_DIR/$BIN_NUMBER
    bin_dir=$output_dir/bins

    mkdir -p $output_dir
    mkdir -p $bin_dir

    bin_length=$((LENGTH / BIN_NUMBER))
    echo "[$(date +"%Y-%m-%d %T")] Simulating $BIN_NUMBER bins with reference length of $LENGTH and bin_length of $bin_length"

    # Simulate reference
    echo "[$(date +"%Y-%m-%d %T")] Simulating genome"
    $BINARY_DIR/mason_genome -l $LENGTH -o $bin_dir/ref.fasta -s $SEED &>/dev/null
    # Evenly distribute it over bins
    echo "[$(date +"%Y-%m-%d %T")] Splitting genome into bins"
    $BINARY_DIR/split_sequence --input $bin_dir/ref.fasta --length $bin_length --parts $BIN_NUMBER
    # We do not need the reference anymore
    rm $bin_dir/ref.fasta

    echo "[$(date +"%Y-%m-%d %T")] Generating $READ_COUNT reads of length $READ_LENGTH with $ERRORS errors"
    read_dir=$output_dir/reads_e$ERRORS\_$READ_LENGTH
    mkdir -p $read_dir
    list_file=$output_dir/${BIN_NUMBER}.filenames
    ln -s $list_file $OUT_DIR/${BIN_NUMBER}.filenames
    seq -f "$output_dir/bins/bin_%0${#BIN_NUMBER}.0f.fasta" 0 1 $((BIN_NUMBER-1)) > $list_file
    $BINARY_DIR/generate_reads \
        --output $read_dir \
        --errors $ERRORS \
        --number_of_reads $READ_COUNT \
        --read_length $READ_LENGTH \
        --number_of_haplotypes $HAPLOTYPE_COUNT \
        $list_file > /dev/null
    find $read_dir -type f -name "*.fastq" -print0 | sort -zV | xargs -0 cat > $read_dir/all
    mv $read_dir/all $read_dir/all.fastq
    # for i in $(seq 0 9); do cat $read_dir/all.fastq >> $read_dir/all_10.fastq; done

    echo ""
done

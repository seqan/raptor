#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -e

READ_LENGTH=100
W=23
K=19
ERRORS=2
HASH=2
SIZES="1g 2g 4g"
THREADS=4
BIN_NUMBER=1024
BINARY_DIR="<path to built binaries>" # containing the raptor binary
INPUT_DIR="<bin path>" # output directory of simulation. the directory that contains the BIN_NUMBER directory
BENCHMARK_DIR="<path>" # directory where results should be stored. E.g., /dev/shm/username; BIN_NUMBER directory will be created.
COPY_INPUT=false # If true, input data will be copied from INPUT_DIR to BENCHMARK_DIR.
EVAL_ENERGY=true # If true, use perf to measure power/energy-pkg/ and power/energy-ram/.

working_directory=$BENCHMARK_DIR/$BIN_NUMBER
mkdir -p $working_directory

if [ "$COPY_INPUT" = true ] ; then
    echo -n "Copying input..."
    mkdir -p $working_directory/bins/
    mkdir -p $working_directory/reads/

    for i in $(seq -f "$INPUT_DIR/$BIN_NUMBER/bins/bin_%0${#BIN_NUMBER}.0f.fasta" 0 1 $((BIN_NUMBER-1)))
    do
        cp $i $working_directory/bins/
    done
    seq -f "$working_directory/bins/bin_%0${#BIN_NUMBER}.0f.fasta" 0 1 $((BIN_NUMBER-1)) > $working_directory/bins.list

    cp $INPUT_DIR/$BIN_NUMBER/reads_e$ERRORS\_$READ_LENGTH/all.fastq $working_directory/reads/
    read_file=$working_directory/reads/all.fastq
    echo "Done."
else
    seq -f "$INPUT_DIR/$BIN_NUMBER/bins/bin_%0${#BIN_NUMBER}.0f.fasta" 0 1 $((BIN_NUMBER-1)) > $working_directory/bins.list
    read_file=$INPUT_DIR/$BIN_NUMBER/reads_e$ERRORS\_$READ_LENGTH/all.fastq
fi

launch_build() {
    if [ "$EVAL_ENERGY" = true ] ; then
        perf stat -o $build_perf -e power/energy-pkg/,power/energy-ram/ "$@"
    else
        "$@"
    fi
}

launch_query() {
    if [ "$EVAL_ENERGY" = true ] ; then
        perf stat -o $query_perf -e power/energy-pkg/,power/energy-ram/ "$@"
    else
        "$@"
    fi
}

for size in $SIZES; do
    ibf_filename=$working_directory/$W\_$K\_$size.ibf # Does not contain HASH
    build_log=$working_directory/$W\_$K\_$size\_build.log
    build_perf=$working_directory/$W\_$K\_$size\_build.perf
    echo "Building IBF with ($W, $K)-minimisers with $HASH hashes and of size $size"
    launch_build    /usr/bin/time -o $build_log -v \
                        $BINARY_DIR/raptor build \
                            --output $ibf_filename \
                            --kmer $K \
                            --window $W \
                            --size $size \
                            --threads $THREADS \
                            --hash $HASH \
                            --input $working_directory/bins.list

    query_log=$working_directory/$W\_$K\_$size\_query.log # Does not contain HASH
    query_perf=$working_directory/$W\_$K\_$size\_query.perf
    query_out=$working_directory/$W\_$K\_$size.out
    echo "Searching IBF for reads of length $READ_LENGTH containing $ERRORS errors"
    launch_query    /usr/bin/time -o $query_log -v \
                            $BINARY_DIR/raptor search \
                                --query $read_file \
                                --index $ibf_filename \
                                --output $query_out \
                                --threads $THREADS \
                                --error $ERRORS \
                                --query_length $READ_LENGTH \
                                --tau 0.9999 \
                                --time

    rm $ibf_filename
done

# Uncomment for basic cleanup, does not delete results
# chmod -R 744 $working_directory/bins
# chmod -R 744 $working_directory/reads
# rm -f $working_directory/bins/*.fasta
# rm -d $working_directory/bins
# rm -f $working_directory/reads/all.fastq
# rm -d $working_directory/reads

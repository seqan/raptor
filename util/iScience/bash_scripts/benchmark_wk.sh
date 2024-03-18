#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -e

READ_LENGTH=250
ERRORS=2
HASH=2
SIZE="4096m"
THREADS=4
BIN_NUMBER=1024
BINARY_DIR="<path to built binaries>" # containing the raptor binary
INPUT_DIR="<bin path>" # output directory of simulation. the directory that contains the BIN_NUMBER directory
BENCHMARK_DIR="<path>" # directory where results should be stored. E.g., /dev/shm/username; BIN_NUMBER directory will be created.
COPY_INPUT=false # If true, input data will be copied from INPUT_DIR to BENCHMARK_DIR.

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

do_task () {
    ibf_filename=$working_directory/$w\_$k\_$SIZE.ibf # Does not contain HASH
    build_log=$working_directory/$w\_$k\_$SIZE\_build.log
    echo "Building IBF with ($w, $k)-minimisers with $HASH hashes and of size $SIZE"
    /usr/bin/time -o $build_log -v \
        $BINARY_DIR/raptor build \
            --output $ibf_filename \
            --kmer $k \
            --window $w \
            --size $SIZE \
            --threads $THREADS \
            --hash $HASH \
            --input $working_directory/bins.list

    query_log=$working_directory/$w\_$k\_$SIZE\_query.log # Does not contain HASH
    query_out=$working_directory/$w\_$k\_$SIZE.out
    echo "Searching IBF for reads of length $READ_LENGTH containing $ERRORS errors"
    /usr/bin/time -o $query_log -v \
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
}

pidlist=""

# 5 jobs in parallel
for w in $(seq 23 2 32 && seq 32 2 80) # w=23,25,27,29,31,32,34,36,38,...,80
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

# Uncomment for basic cleanup, does not delete results
# chmod -R 744 $working_directory/bins
# chmod -R 744 $working_directory/reads
# rm -f $working_directory/bins/*.fasta
# rm -d $working_directory/bins
# rm -f $working_directory/reads/all.fastq
# rm -d $working_directory/reads

#!/usr/bin/env bash
# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2023, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2023, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

set -e

READ_LENGTH=100
W=23
K=19
ERRORS=2
HASH=2
THREADS=32
BIN_NUMBER=65536
BIFROST_BINARY_DIR="/example/build/bifrost/bin" # containing the bifrost binary
INPUT_DIR="/example/big_dataset" # output directory of simulation. the directory that contains the BIN_NUMBER directory
BENCHMARK_DIR="/example/bifrost_bench" # directory where results should be stored. E.g., /dev/shm/username; BIN_NUMBER directory will be created.
COPY_INPUT=true # If true, input data will be copied from INPUT_DIR to BENCHMARK_DIR.

# We need to make bifrost's shared lib available
LD_LIBRARY_PATH=$BIFROST_BINARY_DIR/../lib:$LD_LIBRARY_PATH

working_directory=$BENCHMARK_DIR/$BIN_NUMBER/$FPR
mkdir -p $working_directory

trap 'echo; echo "## [$(date +"%Y-%m-%d %T")] ERROR ##"' ERR

echo "## [$(date +"%Y-%m-%d %T")] Start ##"

if [ "$COPY_INPUT" = true ] ; then
    echo -n "[$(date +"%Y-%m-%d %T")] Copying input..."
    mkdir -p $working_directory/bins/
    mkdir -p $working_directory/reads/

    cp -ur $INPUT_DIR/$BIN_NUMBER/bins $working_directory/

    seq -f "$working_directory/bins/bin_%0${#BIN_NUMBER}.0f.fasta" 0 1 $((BIN_NUMBER-1)) > $working_directory/bins.list

    cp -u $INPUT_DIR/$BIN_NUMBER/reads_e$ERRORS\_$READ_LENGTH/all_10.fastq $working_directory/reads/
    read_file=$working_directory/reads/all_10.fastq
    echo "Done."
else
    seq -f "$INPUT_DIR/$BIN_NUMBER/bins/bin_%0${#BIN_NUMBER}.0f.fasta" 0 1 $((BIN_NUMBER-1)) > $working_directory/bins.list
    read_file=$INPUT_DIR/$BIN_NUMBER/reads_e$ERRORS\_$READ_LENGTH/all_10.fastq
fi

echo -n "[$(date +"%Y-%m-%d %T")] Building Bifrost index..."
build_time=$working_directory/build.time
/usr/bin/time -o $build_time -v \
    $BIFROST_BINARY_DIR/Bifrost build \
        --input-ref-file $working_directory/bins.list \
        --output-file $working_directory/bifrost.index \
        --threads $THREADS \
        --kmer-length $W \
        --min-length $K
echo "Done."

echo -n "[$(date +"%Y-%m-%d %T")] Searching Bifrost index..."
query_time=$working_directory/query.time
/usr/bin/time -o $query_time -v \
    $BIFROST_BINARY_DIR/Bifrost query \
        --input-graph-file $working_directory/bifrost.index.gfa \
        --input-query-file $read_file \
        --output-file $working_directory/query.out \
        --ratio-kmers 0.41 \
        --threads $THREADS
echo "Done."

echo -n "[$(date +"%Y-%m-%d %T")] Cleaning up..."
chmod -R +w $working_directory/bins/
find $working_directory/bins/ -name "*.fasta" -type f -delete
rm -d $working_directory/bins/

chmod -R +w $working_directory/reads/
find $working_directory/reads/ -name "*.fastq" -type f -delete
rm -d $working_directory/reads/
echo "Done."

echo "## [$(date +"%Y-%m-%d %T")] End ##"

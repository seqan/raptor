#!/usr/bin/env bash
# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

set -e

READ_LENGTH=100
W=$3
K=$4
ERRORS=2
HASH=2
PACKSIZES="128 256 512"
FPR="${1}"
THREADS=32
BIN_NUMBER=65536
CHOPPER_BINARY_DIR="/example/build/chopper/bin" # containing the chopper binary
RAPTOR_BINARY_DIR="/example/build/raptor/bin" # containing the raptor binary
INPUT_DIR="/example/big_dataset" # output directory of simulation. the directory that contains the BIN_NUMBER directory
BENCHMARK_DIR="/example/raptor_bench" # directory where results should be stored. E.g., /dev/shm/username; BIN_NUMBER directory will be created.
COPY_INPUT=true # If true, input data will be copied from INPUT_DIR to BENCHMARK_DIR.

working_directory=$BENCHMARK_DIR/$BIN_NUMBER/$FPR
mkdir -p $working_directory
mkdir -p $working_directory/hll

trap 'echo; echo "## [$(date +"%Y-%m-%d %T")] ERROR ##"' ERR

echo "## [$(date +"%Y-%m-%d %T")] Start: FPR=$FPR, ($W, $K)-minimizer, IBF=$2 ##"

if [ "$COPY_INPUT" = true ] ; then
    echo -n "[$(date +"%Y-%m-%d %T")] Copying input..."
    mkdir -p $working_directory/bins/
    mkdir -p $working_directory/reads/

    cp -ur $INPUT_DIR/$BIN_NUMBER/bins $working_directory/

    seq -f "$working_directory/bins/bin_%0${#BIN_NUMBER}g.fasta" 0 1 $((BIN_NUMBER-1)) > $working_directory/bins.list

    cp -u $INPUT_DIR/$BIN_NUMBER/reads_e$ERRORS\_$READ_LENGTH/all_10.fastq $working_directory/reads/
    read_file=$working_directory/reads/all_10.fastq
    echo "Done."
else
    seq -f "$INPUT_DIR/$BIN_NUMBER/bins/bin_%0${#BIN_NUMBER}g.fasta" 0 1 $((BIN_NUMBER-1)) > $working_directory/bins.list
    read_file=$INPUT_DIR/$BIN_NUMBER/reads_e$ERRORS\_$READ_LENGTH/all_10.fastq
fi

echo -n "[$(date +"%Y-%m-%d %T")] Counting k-mers..."
count_file=$working_directory/$W\_$K\_all_bins.counts
count_time=$working_directory/$W\_$K\_count.time
/usr/bin/time -o $count_time -v \
    $CHOPPER_BINARY_DIR/chopper count \
        --exclusively-hlls \
        -t $THREADS \
        -c 1 \
        --hll-dir $working_directory/hll \
        --data_file $working_directory/bins.list \
        --outfile $count_file
echo "Done."

best_pack_file=$working_directory/$W\_$K\_best_pack.pack
best_pack_time=$working_directory/$W\_$K\_best_pack.time
echo "[$(date +"%Y-%m-%d %T")] Determining best t_max..."
/usr/bin/time -o $best_pack_time -v \
    $CHOPPER_BINARY_DIR/chopper pack \
        -f $count_file \
        -b 1024 \
        --hll-dir $working_directory/hll \
        --estimate-union \
        --determine-num-bins \
        --num-threads $THREADS \
        -o $best_pack_file \
        &>$best_pack_file.txt
echo "Done."

for pack in $PACKSIZES; do
    pack_file=$working_directory/$W\_$K\_$pack.pack
    pack_time=$working_directory/$W\_$K\_$pack\_pack.time
    echo -n "[$(date +"%Y-%m-%d %T")] Computing layout for index with ($W, $K)-minimisers, $HASH hashes, t_max=$pack..."
    /usr/bin/time -o $pack_time -v \
        $CHOPPER_BINARY_DIR/chopper pack \
            -f $count_file \
            -b $pack \
            --hll-dir $working_directory/hll \
            --estimate-union \
            --num-threads $THREADS \
            -o $pack_file \
            1>/dev/null
    echo "Done."

    index_filename=$working_directory/$W\_$K\_$pack.index # Does not contain HASH
    build_time=$working_directory/$W\_$K\_$pack\_build.time
    echo -n "[$(date +"%Y-%m-%d %T")] Building index with ($W, $K)-minimisers, $HASH hashes, t_max=$pack..."
    /usr/bin/time -o $build_time -v \
        $RAPTOR_BINARY_DIR/raptor build \
            --hibf \
            --fpr $FPR \
            --output $index_filename \
            --kmer $K \
            --window $W \
            --threads $THREADS \
            --hash $HASH \
            $pack_file
    echo "Done."

    query_out=$working_directory/$W\_$K\_$pack.out # Does not contain HASH
    query_time=$working_directory/$W\_$K\_$pack\_query.time
    echo -n "[$(date +"%Y-%m-%d %T")] Searching index for reads of length $READ_LENGTH containing $ERRORS errors..."
    /usr/bin/time -o $query_time -v \
        $RAPTOR_BINARY_DIR/raptor search \
            --hibf \
            --query $read_file \
            --index $index_filename \
            --output $query_out \
            --threads $THREADS \
            --error $ERRORS \
            --pattern $READ_LENGTH \
            --tau 0.95 \
            --time
    echo "Done."
done

### VANILLA RAPTOR START

index_filename=$working_directory/$W\_$K\_IBF.index # Does not contain HASH
build_time=$working_directory/$W\_$K\_IBF\_build.time
echo -n "[$(date +"%Y-%m-%d %T")] Building index with ($W, $K)-minimisers, $HASH hashes, size=${2}..."
/usr/bin/time -o $build_time -v \
    $RAPTOR_BINARY_DIR/raptor build \
        --output $index_filename \
        --kmer $K \
        --window $W \
        --size "${2}" \
        --threads $THREADS \
        --hash $HASH \
        $working_directory/bins.list

echo "Done."

query_out=$working_directory/$W\_$K\_IBF.out # Does not contain HASH
query_time=$working_directory/$W\_$K\_IBF\_query.time
echo -n "[$(date +"%Y-%m-%d %T")] Searching index for reads of length $READ_LENGTH containing $ERRORS errors..."
/usr/bin/time -o $query_time -v \
    $RAPTOR_BINARY_DIR/raptor search \
        --query $read_file \
        --index $index_filename \
        --output $query_out \
        --threads $THREADS \
        --error $ERRORS \
        --pattern $READ_LENGTH \
        --tau 0.95 \
        --time
echo "Done."

### VANILLA RAPTOR END

echo -n "[$(date +"%Y-%m-%d %T")] Cleaning up..."
chmod -R +w $working_directory/bins/
find $working_directory/bins/ -name "*.fasta" -type f -delete
rm -d $working_directory/bins/

chmod -R +w $working_directory/reads/
find $working_directory/reads/ -name "*.fastq" -type f -delete
rm -d $working_directory/reads/
echo "Done."

echo "## [$(date +"%Y-%m-%d %T")] End ##"

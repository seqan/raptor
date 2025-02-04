#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -e

ERRORS=2
HASH=2
THREADS=32
BINARY_DIR="<path to built binaries>" # containing the DREAM-Yara binaries
BIN_DIR="<bin path>" # directory containing the binned reference data set
BENCHMARK_DIR="<path>" # directory where the benchmarks should be run. Input data will be copied here if COPY_INPUT is true. e.g. /dev/shm/username; BIN_NUMBER directory will be created.
COPY_INPUT=false # copy the input (bins, reads) to the BENCHMARK_DIR. useful when the BENCHMARK_DIR is set to /dev/shm/

do_task () {
    working_directory=$BENCHMARK_DIR/$BIN_NUMBER
    mkdir -p $working_directory/fm_indices/

    identifier=$K\_$SIZE\G

    if [ "$COPY_INPUT" = true ] ; then
        echo "[$(date +%H:%M)] Copying input"
        mkdir -p $working_directory/bins/
        mkdir -p $working_directory/reads/
        for i in $(seq -f "$BIN_DIR/$BIN_NUMBER/bins/bin_%0${#BIN_NUMBER}.0f.fasta.gz" 0 1 $((BIN_NUMBER-1))); do
            cp $i $working_directory/bins/
        done
        cp $BIN_DIR/$BIN_NUMBER/reads_e$ERRORS/$READ_LENGTH/all.fastq $working_directory/reads/
        bin_list=$(seq -f "$working_directory/bins/bin_%0${#BIN_NUMBER}.0f.fasta.gz" 0 1 $((BIN_NUMBER-1)))
        read_file=$working_directory/reads/all.fastq
        COPY_INPUT=false
    else
        bin_list=$(seq -f "$BIN_DIR/$BIN_NUMBER/bins/bin_%0${#BIN_NUMBER}.0f.fasta.gz" 0 1 $((BIN_NUMBER-1)))
        read_file=$BIN_DIR/$BIN_NUMBER/reads_e$ERRORS/$READ_LENGTH/all.fastq
    fi

    ####################################################################################################################
    ################################################### Build IBF ###################################################
    ####################################################################################################################
    ibf_filename=$working_directory/$identifier.filter
    build_ibf_time=$working_directory/$identifier\_build_ibf.time
    build_ibf_log=$working_directory/$identifier\_build_ibf.log

    if [ "$BUILD_IBF" = true ]; then
        echo "[$(date +%H:%M)] Building IBF with $K-mers and a size of $SIZE GiB"
        /usr/bin/time -o $build_ibf_time -v \
            $BINARY_DIR/dream_yara_build_filter \
                --output-file $ibf_filename \
                --kmer-size $K \
                --bloom-size $SIZE \
                --threads $THREADS \
                --num-hash $HASH \
                --verbose \
                --version-check 0 \
                $bin_list \
                &> $build_ibf_log
    else
        echo "[$(date +%H:%M)] Skipping building of IBF with $K-mers and a size of $SIZE GiB"
    fi

    ####################################################################################################################
    ################################################# Build FM-indices #################################################
    ####################################################################################################################
    build_fm_indices_time=$working_directory/$identifier\_build_fm.time
    build_fm_indices_log=$working_directory/$identifier\_build_fm.log

    if [ "$BUILD_INDEX" = true ] ; then
        echo "[$(date +%H:%M)] Building FM-indices for $BIN_NUMBER bins"
        /usr/bin/time -o $build_fm_indices_time -v \
            $BINARY_DIR/dream_yara_indexer \
                --output-prefix $working_directory/fm_indices/ \
                --threads $THREADS \
                --verbose \
                --version-check 0 \
                $bin_list \
                &> $build_fm_indices_log

        BUILD_INDEX=false
    else
        echo "[$(date +%H:%M)] Skipping building of FM-Indices for $BIN_NUMBER bins"
    fi

    ####################################################################################################################
    ################################################# Search for query #################################################
    ####################################################################################################################
    mapper_time=$working_directory/$identifier\_mapper\_$READ_LENGTH.time
    mapper_log=$working_directory/$identifier\_mapper\_$READ_LENGTH.log
    mapper_out=$working_directory/$identifier\_$READ_LENGTH.sam

    echo "[$(date +%H:%M)] Mapping reads of length $READ_LENGTH containing $ERRORS + 1 errors"
    until /usr/bin/time -o $mapper_time -v \
        $BINARY_DIR/dream_yara_mapper \
            --bloom-filter $ibf_filename \
            --output-file $mapper_out \
            --threads $THREADS \
            --error-rate 0$(bc -l <<< "($ERRORS+1)/$READ_LENGTH") \
            --verbose \
            --version-check 0 \
            $working_directory/fm_indices/ \
            $read_file \
            &> $mapper_log
    do
        sleep 5
    done

    # rm $ibf_filename
}

for BIN_NUMBER in 1024 64; do
    BUILD_INDEX=false
    for SIZE in 16 8; do
        BUILD_IBF=true
        for READ_LENGTH in 100 250; do
            K=19
            do_task

            K=23
            do_task

            K=31
            do_task

            BUILD_IBF=false
        done
    done
done

# Uncomment for basic cleanup, does not delete results
# chmod -R 744 $working_directory/bins
# chmod -R 744 $working_directory/reads
# rm -f $working_directory/bins/*.fasta
# rm -d $working_directory/bins
# rm -f $working_directory/reads/all.fastq
# rm -d $working_directory/reads

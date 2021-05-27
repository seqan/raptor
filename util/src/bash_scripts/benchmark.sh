#!/bin/bash
set -e

READ_LENGTH=250
ERRORS=2
HASH=2
SIZE="4096m"
THREADS=4
BIN_NUMBER=1024
BINARY_DIR="<path to built binaries>" # containing the raptor binary
BIN_DIR="<bin path>" # output directory of simulation. the directory that contains the BIN_NUMBER directory
BENCHMARK_DIR="<path>" # directory where the benchmarks should be run. Input data will be copied here. e.g. /dev/shm/username; BIN_NUMBER directory will be created.

working_directory=$BENCHMARK_DIR/$BIN_NUMBER
mkdir -p $working_directory/bins/
mkdir -p $working_directory/reads/

for i in $(seq -f "$BIN_DIR/$BIN_NUMBER/bins/bin_%0${#BIN_NUMBER}g.fasta" 0 1 $((BIN_NUMBER-1)))
do
    cp $i $working_directory/bins/
done
cp $BIN_DIR/$BIN_NUMBER/reads_e$ERRORS\_$READ_LENGTH/all.fastq $working_directory/reads/

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
            $(seq -f "$working_directory/bins/bin_%0${#BIN_NUMBER}g.fasta" 0 1 $((BIN_NUMBER-1)))

    query_log=$working_directory/$w\_$k\_$SIZE\_query.log # Does not contain HASH
    query_out=$working_directory/$w\_$k\_$SIZE.out
    echo "Searching IBF for reads of length $READ_LENGTH containing $ERRORS errors"
    /usr/bin/time -o $query_log -v \
        $BINARY_DIR/raptor search \
            --query $working_directory/reads/all.fastq \
            --index $ibf_filename \
            --output $query_out \
            --kmer $k \
            --window $w \
            --threads $THREADS \
            --error $ERRORS \
            --pattern $READ_LENGTH \
            --tau 0.9999 \
            --time

    rm $ibf_filename
}

pidlist=""

for w in $(seq 23 2 32 && seq 32 2 80)
do
    for k in 16 17 18 19 20
    do
        if [[ $w != 23 ]] || [[ $k != 20 ]]; then # Segfault for READ_LENGTH = 150
            do_task & pidlist="$pidlist $!"
        fi
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

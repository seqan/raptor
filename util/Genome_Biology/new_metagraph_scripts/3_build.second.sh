#!/usr/bin/env bash
set -Eeuo pipefail

source /project/archive-index-data/seiler/refseq_metagraph/variables.sh

/usr/bin/time -v -o ${TIME_DIR}/build.second.time \
    /project/archive-index-data/software/bin/metagraph build \
        --state fast \
        --mode canonical \
        --parallel ${NUM_THREADS} \
        -k ${KMER_SIZE} \
        --mem-cap-gb 10 \
        -o ${INDEX_DIR}/second \
        ${INDEX_DIR}/transformed.fasta.gz \
    &> ${LOG_DIR}/build.second.log

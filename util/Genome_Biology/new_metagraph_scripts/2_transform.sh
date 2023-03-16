#!/usr/bin/env bash
set -Eeuo pipefail

source /project/archive-index-data/seiler/refseq_metagraph/variables.sh

/usr/bin/time -v -o ${TIME_DIR}/transform.time \
    /project/archive-index-data/software/bin/metagraph transform \
        --to-fasta \
        --primary-kmers \
        -p ${NUM_THREADS} \
        -o ${INDEX_DIR}/transformed \
        ${INDEX_DIR}/first.dbg \
    &> ${LOG_DIR}/transform.log

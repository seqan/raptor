#!/usr/bin/env bash
set -Eeuo pipefail

source /project/archive-index-data/seiler/refseq_metagraph/variables.sh

cat ${KMC_FILENAMES} | \
    /usr/bin/time -v -o ${TIME_DIR}/build.first.time \
        /project/archive-index-data/software/bin/metagraph build \
            --state fast \
            --mode canonical \
            --parallel ${NUM_THREADS} \
            -k ${KMER_SIZE} \
            --mem-cap-gb 10 \
            -o ${INDEX_DIR}/first \
        &> ${LOG_DIR}/build.first.log

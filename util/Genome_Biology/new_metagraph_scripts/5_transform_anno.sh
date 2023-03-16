#!/usr/bin/env bash
set -Eeuo pipefail

source /project/archive-index-data/seiler/refseq_metagraph/variables.sh

find ${ANNO_DIR} -name *.column.annodbg | \
    /usr/bin/time -v -o ${TIME_DIR}/transform_anno.time \
        /project/archive-index-data/software/bin/metagraph transform_anno \
            -p ${NUM_THREADS} \
            --anno-type brwt \
            --greedy \
            --fast \
            -o ${INDEX_DIR}/transformed
        &> ${LOG_DIR}/transform_anno.log

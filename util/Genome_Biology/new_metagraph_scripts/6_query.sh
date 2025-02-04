#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -Eeuo pipefail

source /project/archive-index-data/seiler/refseq_metagraph/variables.sh

QUERY_LIST="10M 5M 1M 1K"

for QUERY in ${QUERY_LIST}; do
    QUERY_FILE="QUERY_${QUERY}"

    /usr/bin/time -v -o ${TIME_DIR}/query.${QUERY}.time \
        /project/archive-index-data/software/bin/metagraph query \
            -i ${INDEX_DIR}/second.dbg \
            -a ${INDEX_DIR}/transformed.brwt.annodbg \
            --discovery-fraction ${QUERY_THRESHOLD} \
            --parallel ${NUM_THREADS} \
            --fast \
            ${!QUERY_FILE} \
            > ${WORK_DIR}/query.${QUERY}.result \
            2> ${LOG_DIR}/query.${QUERY}.log
done

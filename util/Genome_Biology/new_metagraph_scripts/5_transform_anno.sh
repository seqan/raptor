#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

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

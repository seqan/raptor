#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -Eeuo pipefail

source /project/archive-index-data/seiler/refseq_metagraph/variables.sh

cat ${KMC_FILENAMES} | \
    /usr/bin/time -v -o ${TIME_DIR}/annotate.time \
        /project/archive-index-data/software/bin/metagraph annotate \
            -i ${INDEX_DIR}/second.dbg \
            --anno-filename \
            --separately \
            -p ${NUM_THREADS} \
            --threads-each 1 \
            -o ${ANNO_DIR} \
        &> ${LOG_DIR}/annotate.log

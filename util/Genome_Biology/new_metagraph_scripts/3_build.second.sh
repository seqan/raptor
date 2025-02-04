#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

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

#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

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

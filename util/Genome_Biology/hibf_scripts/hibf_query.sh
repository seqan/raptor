#!/bin/bash

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

source /project/archive-index-data/smehringer/benchmark.variables

RAPTOR_DIR="${WORKDIR}/raptor_bench_new"

/project/archive-index-data/software/bin/raptor search --index ${RAPTOR_DIR}/raptor.index \
                      --query ${QUERY_FILE} \
                      --output ${RAPTOR_DIR}/raptor.result \
                      --threshold ${QUERY_THRESHOLD} \
                      --threads ${NUM_THREADS}

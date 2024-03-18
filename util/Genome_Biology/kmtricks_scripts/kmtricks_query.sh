# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

source /project/archive-index-data/smehringer/benchmark.variables

KM_DIR="${WORKDIR}/kmtricks_bench"

# Query

#/project/archive-index-data/software/bin/kmtricks query --run-dir ${KM_DIR}/kmtricks_index --query ${QUERY_FILE_FASTA} --threshold ${QUERY_THRESHOLD} --no-detail --threads ${NUM_THREADS} > ${KM_DIR}/kmtricks.result
/project/archive-index-data/software/bin/kmtricks query --run-dir ${KM_DIR}/kmtricks_index --query ${QUERY_FILE_FASTA} --threshold 0.6 --no-detail --threads ${NUM_THREADS} > ${KM_DIR}/kmtricks.result

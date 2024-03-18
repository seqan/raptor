# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

source /project/archive-index-data/smehringer/benchmark.variables

COBS_DIR="${WORKDIR}/COBS_compact_bench"

# Query

/project/archive-index-data/software/bin/cobs query --index ${COBS_DIR}/index.cobs_compact --file ${QUERY_FILE_FASTA} --threads ${NUM_THREADS} --threshold ${QUERY_THRESHOLD} > ${COBS_DIR}/COBS.result

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

source /project/archive-index-data/smehringer/benchmark.variables

BIF_DIR="${WORKDIR}/bifrost_bench"

/project/archive-index-data/software/bin/Bifrost query \
        --input-graph-file ${BIF_DIR}/bifrost.index.gfa.gz \
        --input-query-file ${QUERY_FILE} \
        --output-file ${BIF_DIR}/bifrost.result \
        --ratio-kmers ${QUERY_THRESHOLD} \
        --input-color-file ${BIF_DIR}/bifrost.index.color.bfg \
        --threads ${NUM_THREADS} 

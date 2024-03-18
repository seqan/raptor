# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

source /project/archive-index-data/smehringer/benchmark.variables

BIF_DIR="${WORKDIR}/bifrost_bench"

/project/archive-index-data/software/bin/Bifrost build --input-ref-file ${DATA_FILENAMES} --output-file ${BIF_DIR}/bifrost.index  --kmer-length 31 --threads ${NUM_THREADS} --colors

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

source /project/archive-index-data/smehringer/benchmark.variables

RAPTOR_DIR="${WORKDIR}/ibf_bench"

/project/archive-index-data/software/bin/raptor build --threads ${NUM_THREADS} --kmer ${KMER_SIZE} --fpr 0.05 --hash 2 ${DATA_FILENAMES} --output ${RAPTOR_DIR}/raptor.index

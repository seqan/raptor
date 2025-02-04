# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

source /project/archive-index-data/smehringer/benchmark.variables

META_DIR="${WORKDIR}/metagraph_bench"

# Metagraph build

cat ${DATA_FILENAMES} | /project/archive-index-data/software/bin/metagraph build --state fast --mode canonical --parallel ${NUM_THREADS} -k ${KMER_SIZE} --mem-cap-gb 10 -o ${META_DIR}/metagraph_index/graph 2>&1 | tee ${META_DIR}/build_log.txt


# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

source /project/archive-index-data/smehringer/benchmark.variables

META_DIR="${WORKDIR}/metagraph_bench"

# Annotate

cat ${DATA_FILENAMES} | /project/archive-index-data/software/bin/metagraph annotate -v --parallel ${NUM_THREADS} -i ${META_DIR}/metagraph_index/graph.dbg --mem-cap-gb 10 --anno-filename -o ${META_DIR}/metagraph_index/annotation


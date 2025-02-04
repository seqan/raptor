# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

source /project/archive-index-data/smehringer/benchmark.variables

KM_DIR="${WORKDIR}/kmtricks_bench"

echo "KMTRICKS INDEX"
/project/archive-index-data/software/bin/kmtricks index --run-dir ${KM_DIR}/kmtricks_index --howde --threads ${NUM_THREADS}



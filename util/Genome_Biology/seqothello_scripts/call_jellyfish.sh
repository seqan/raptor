#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -Eeuo pipefail

source /project/archive-index-data/smehringer/benchmark.variables

OTH_DIR="${WORKDIR}/SeqOthello_bench"

INPUT="$1"
FILENAME=$(basename ${INPUT})

/project/archive-index-data/software/bin/jellyfish count -s 100M -m ${KMER_SIZE} -C -o ${OTH_DIR}/kmers/${FILENAME}.jf <(zcat ${INPUT})
/project/archive-index-data/software/bin/jellyfish dump -t -L 1 -c -o ${OTH_DIR}/kmers/${FILENAME}.kmer ${OTH_DIR}/kmers/${FILENAME}.jf

rm ${OTH_DIR}/kmers/${FILENAME}.jf

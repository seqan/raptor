# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

source /project/archive-index-data/smehringer/benchmark.variables

OTH_DIR="${WORKDIR}/SeqOthello_bench"

# Query

/project/archive-index-data/software/bin/seqothello-query --map-folder=${OTH_DIR}/SeqOthello_index/ --transcript=${QUERY_FILE_FASTA} --qthread=${NUM_THREADS} --output=${OTH_DIR}/SeqOthello.result

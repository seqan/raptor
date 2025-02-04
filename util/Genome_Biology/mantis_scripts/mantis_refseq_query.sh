#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

source /project/archive-index-data/smehringer/benchmark.variables

LINES=$(grep @ ${QUERY_FILE} | wc -l )
LINES2=$(cat /project/archive-index-data/smehringer/mantis_bench/mantis.query | wc -l)

if [[ "$LINES" != "$LINES2" ]]
then
	echo "Number of queries doesn't match. You might have to execute copy_input.sh"
	exit 1
fi


/project/archive-index-data/software/bin/mantis-master query \
    -p /project/archive-index-data/smehringer/mantis_bench/mantis/ \
    -k ${KMER_SIZE} \
    -o /project/archive-index-data/smehringer/mantis_bench/mantis.result \
    /project/archive-index-data/smehringer/mantis_bench/mantis.query

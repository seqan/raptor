#!/usr/bin/env bash
# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2021, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2021, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

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

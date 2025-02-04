#!/bin/bash

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -eEuo pipefail

source /path/to/benchmark.variables

LOG_FILE="]${WORKDIR}/run_all_benchmarks.log"
TIMINGS="${WORKDIR}/all_benchmarks.timings"
MEMORY="${WORKDIR}/all_benchmarks.mem"
INDEX_SIZES="${WORKDIR}/all_benchmarks.index_size"

scripts=( \
# Build
${WORKDIR}/PAC_bench/PAC_build_me.sh \
${WORKDIR}/metagraph_bench/metagraph_build.sh \
${WORKDIR}/metagraph_bench/metagraph_annotate.sh \
${WORKDIR}/COBS_compact_bench/COBS_build.sh \
${WORKDIR}/SeqOthello_bench/SeqOthello_prep.sh \
${WORKDIR}/SeqOthello_bench/SeqOthello_build.sh \
${WORKDIR}/mantis_bench/mantis_build.sh \
${WORKDIR}/ibf_bench/raptor_build.sh \
${WORKDIR}/bifrost_bench/bifrost_build.sh \
${WORKDIR}/hibf_bench/hibf_build.sh
# Query
${WORKDIR}/COBS_compact_bench/COBS_query.sh \
${WORKDIR}/COBS_compact_bench/COBS_query.sh \
${WORKDIR}/COBS_compact_bench/COBS_query.sh \
${WORKDIR}/metagraph_bench/metagraph_query.sh \
${WORKDIR}/metagraph_bench/metagraph_query.sh \
${WORKDIR}/metagraph_bench/metagraph_query.sh \
${WORKDIR}/SeqOthello_bench/SeqOthello_query.sh \
${WORKDIR}/SeqOthello_bench/SeqOthello_query.sh \
${WORKDIR}/SeqOthello_bench/SeqOthello_query.sh \
${WORKDIR}/hibf_bench/hibf_query.sh \
${WORKDIR}/hibf_bench/hibf_query.sh \
${WORKDIR}/hibf_bench/hibf_query.sh \
${WORKDIR}/mantis_bench/mantis_query.sh \
${WORKDIR}/mantis_bench/mantis_query.sh \
${WORKDIR}/mantis_bench/mantis_query.sh \
${WORKDIR}/ibf_bench/raptor_query.sh \
${WORKDIR}/ibf_bench/raptor_query.sh \
${WORKDIR}/ibf_bench/raptor_query.sh \
${WORKDIR}/bifrost_bench/bifrost_query.sh \
${WORKDIR}/bifrost_bench/bifrost_query.sh \
${WORKDIR}/bifrost_bench/bifrost_query.sh \
 )

# Run Prep

for SH_FILE in "${scripts[@]}"
do
        LOCAL_LOG_FILE="/tmp/smehringer/local.log"
        echo "RUNNING ${SH_FILE}" | tee -a ${LOG_FILE}

        /usr/bin/time -o ${LOCAL_LOG_FILE} -v ${SH_FILE} &>> ${LOG_FILE}
        #/usr/bin/time -o ${LOCAL_LOG_FILE} -v ./echo &>> ${LOG_FILE}

        echo -en "${SH_FILE}\t" >> ${TIMINGS}
        echo -en "${SH_FILE}\t" >> ${MEMORY}
        grep Elapsed ${LOCAL_LOG_FILE} | cut -d ' ' -f 8 >> ${TIMINGS}
        grep Max ${LOCAL_LOG_FILE} | cut -d ' ' -f 6 | numfmt --to iec --format "%8.2f" --from-unit Ki --round nearest >> ${MEMORY}

        cat ${LOCAL_LOG_FILE} >> ${LOG_FILE}
done

# collect index sizes
#for INDEX in $(find -maxdepth 2 -path "*_bench/*index*"); do du -sh ${INDEX} >> ${INDEX_SIZES}; done

echo "DONE."

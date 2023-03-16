#!/usr/bin/env bash
set -Eeuo pipefail

source /project/archive-index-data/seiler/refseq_metagraph/variables.sh

mkdir -p ${WORK_DIR}/kmc_tmp

# Optionally: kmc spawns 2 threads
# NUM_THREADS=$((NUM_THREADS/2))

/usr/bin/time -v -o ${TIME_DIR}/kmc.time \
    parallel \
        --line-buffer \
        -j${NUM_THREADS} \
        -a "${DATA_FILENAMES}" \
        $(dirname "$0")/kmc.helper.sh \
    &> ${LOG_DIR}/kmc.log

find ${KMC_DIR} -name "*.kmc_suf" > ${KMC_FILENAMES}

rm -d ${WORK_DIR}/kmc_tmp

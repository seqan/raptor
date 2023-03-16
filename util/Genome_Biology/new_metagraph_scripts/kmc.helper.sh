#!/usr/bin/env bash
set -Eeuo pipefail

source /project/archive-index-data/seiler/refseq_metagraph/variables.sh

INPUT="$1"
OUTPUT="${KMC_DIR}/$(basename ${INPUT} .fna.gz)"

/project/archive-index-data/software/bin/kmc \
    -fm \
    -t1 \
    -r \
    -k${KMER_SIZE} \
    -ci1 \
    -hp \
    ${INPUT} \
    ${OUTPUT} \
    ${WORK_DIR}/kmc_tmp \
    1> /dev/null

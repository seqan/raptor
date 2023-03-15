#!/usr/bin/env bash
set -Eeuo pipefail

source /project/archive-index-data/smehringer/benchmark.variables

OTH_DIR="${WORKDIR}/SeqOthello_bench"

INPUT="$1"
FILENAME=$(basename ${INPUT})

/project/archive-index-data/software/bin/jellyfish count -s 100M -m ${KMER_SIZE} -C -o ${OTH_DIR}/kmers/${FILENAME}.jf <(zcat ${INPUT})
/project/archive-index-data/software/bin/jellyfish dump -t -L 1 -c -o ${OTH_DIR}/kmers/${FILENAME}.kmer ${OTH_DIR}/kmers/${FILENAME}.jf

rm ${OTH_DIR}/kmers/${FILENAME}.jf

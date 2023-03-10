#!/usr/bin/env bash
set -Eeuo pipefail

source /project/archive-index-data/smehringer/benchmark.variables


experiment="$1"
#FILENAME=$(basename ${INPUT})



HOW_DIR="${WORKDIR}/howdesbt_bench"

gzip -dc ${experiment} > ${HOW_DIR}/$(basename ${experiment} .gz) 
/project/archive-index-data/software/bin/howdesbt makebf K=${KMER_SIZE} --bits=125488049 ${HOW_DIR}/$(basename ${experiment} .gz) --out=${HOW_DIR}/bfs/$(basename ${experiment} fna.gz).bf
rm ${HOW_DIR}/$(basename ${experiment} .gz)


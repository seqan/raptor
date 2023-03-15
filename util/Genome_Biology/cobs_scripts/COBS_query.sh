source /project/archive-index-data/smehringer/benchmark.variables

COBS_DIR="${WORKDIR}/COBS_compact_bench"

# Query

/project/archive-index-data/software/bin/cobs query --index ${COBS_DIR}/index.cobs_compact --file ${QUERY_FILE_FASTA} --threads ${NUM_THREADS} --threshold ${QUERY_THRESHOLD} > ${COBS_DIR}/COBS.result

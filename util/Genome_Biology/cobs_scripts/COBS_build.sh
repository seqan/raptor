source /project/archive-index-data/smehringer/benchmark.variables

COBS_DIR="${WORKDIR}/COBS_compact_bench"

cp ${DATA_FILENAMES} ${COBS_DIR}/data.list

# BUild
/project/archive-index-data/software/bin/cobs compact-construct --false-positive-rate ${FPR}    \
                                                                --num-hashes 2                  \
                                                                -k ${KMER_SIZE}                 \
                                                                --threads ${NUM_THREADS}        \
                                                                --file-type list                \
                                                                ${COBS_DIR}/data.list           \
                                                                ${COBS_DIR}/index.cobs_compact

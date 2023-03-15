source /project/archive-index-data/smehringer/benchmark.variables

KM_DIR="${WORKDIR}/kmtricks_bench"

# Query

#/project/archive-index-data/software/bin/kmtricks query --run-dir ${KM_DIR}/kmtricks_index --query ${QUERY_FILE_FASTA} --threshold ${QUERY_THRESHOLD} --no-detail --threads ${NUM_THREADS} > ${KM_DIR}/kmtricks.result
/project/archive-index-data/software/bin/kmtricks query --run-dir ${KM_DIR}/kmtricks_index --query ${QUERY_FILE_FASTA} --threshold 0.6 --no-detail --threads ${NUM_THREADS} > ${KM_DIR}/kmtricks.result

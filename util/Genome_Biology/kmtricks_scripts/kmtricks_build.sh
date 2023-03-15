source /project/archive-index-data/smehringer/benchmark.variables

KM_DIR="${WORKDIR}/kmtricks_bench"

echo "KMTRICKS INDEX"
/project/archive-index-data/software/bin/kmtricks index --run-dir ${KM_DIR}/kmtricks_index --howde --threads ${NUM_THREADS}



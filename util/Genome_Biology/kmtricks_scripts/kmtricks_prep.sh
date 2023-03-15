source /project/archive-index-data/smehringer/benchmark.variables

KM_DIR="${WORKDIR}/kmtricks_bench"

# Collect data 

grep "^#" ${WORKDIR}/raptor_bench/the_one_and_only.truth | grep -v QUERY_NAME | tr '#' 'D' | sed -e 's/\t/: /g' > ${KM_DIR}/data.fof

# Build
echo "KMTRICKS PIPELINE"
/project/archive-index-data/software/bin/kmtricks pipeline --file ${KM_DIR}/data.fof             \
                                                           --run-dir ${KM_DIR}/kmtricks_index    \
                                                           --kmer-size ${KMER_SIZE}              \
                                                           --mode hash:bft:bin         \
                                                           --hard-min 1                \
                                                           --soft-min 1                \
                                                           --share-min 1               \
                                                           --bloom-size 125488049      \
                                                           --bf-format howdesbt        \
                                                           --threads 16                \
                                                           --cpr                       \
                                                           --skip-merge


source /project/archive-index-data/smehringer/benchmark.variables 

RAPTOR_DIR="${WORKDIR}/ibf_bench"

/project/archive-index-data/software/bin/raptor build --threads ${NUM_THREADS} --kmer ${KMER_SIZE} --fpr 0.05 --hash 2 ${DATA_FILENAMES} --output ${RAPTOR_DIR}/raptor.index

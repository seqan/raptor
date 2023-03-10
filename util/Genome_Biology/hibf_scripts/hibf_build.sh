source /project/archive-index-data/smehringer/benchmark.variables 

RAPTOR_DIR="${WORKDIR}/raptor_bench_new"

/project/archive-index-data/software/bin/raptor layout --tmax 192 --input-file ${DATA_FILENAMES} --kmer-size ${KMER_SIZE} --output-filename ${RAPTOR_DIR}/raptor.layout --false-positive-rate 0.05 --threads ${NUM_THREADS} --rearrange-user-bins --num-hash-functions 2
/project/archive-index-data/software/bin/raptor build --threads ${NUM_THREADS} --kmer ${KMER_SIZE} --fpr 0.05 --hash 2 ${RAPTOR_DIR}/raptor.layout --output ${RAPTOR_DIR}/raptor.index

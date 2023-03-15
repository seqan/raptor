source /project/archive-index-data/smehringer/benchmark.variables

OTH_DIR="${WORKDIR}/SeqOthello_bench"

# Query

/project/archive-index-data/software/bin/seqothello-query --map-folder=${OTH_DIR}/SeqOthello_index/ --transcript=${QUERY_FILE_FASTA} --qthread=${NUM_THREADS} --output=${OTH_DIR}/SeqOthello.result

source /project/archive-index-data/smehringer/benchmark.variables

META_DIR="${WORKDIR}/metagraph_bench"

# Query MetaGraph

/project/archive-index-data/software/bin/metagraph query -v -i ${META_DIR}/metagraph_index/graph.dbg -a ${META_DIR}/metagraph_index/annotation.column.annodbg --discovery-fraction ${QUERY_THRESHOLD} --parallel ${NUM_THREADS} --fast ${QUERY_FILE} > ${META_DIR}/metagraph.result 

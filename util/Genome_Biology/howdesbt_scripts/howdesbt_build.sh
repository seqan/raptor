# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

source /project/archive-index-data/smehringer/benchmark.variables

HOW_DIR="${WORKDIR}/howdesbt_bench"

# Following https://github.com/medvedevgroup/HowDeSBT/tree/master/tutorial

# Bloom filter size chosen for 5% FPR, ont he biggest RefSeq Count of 15880273 32-kmers
# via scipt:
# python2 simple_bf_size_estimate.py 15880273 5%
# #numItems	bfFP	numHashes	numBits
# 15880273	0.050000	2	125488049


# Create Bloom Filters
#mkdir -p ${HOW_DIR}/bfs

#parallel --line-buffer -j ${NUM_THREADS} -a ${DATA_FILENAMES} ${HOW_DIR}/call_makebf.sh

# CLuster leaved
# bits: Use 10% of the original bit size as
find ${HOW_DIR}/bfs/ -name *.bf > ${HOW_DIR}/leafnames
/project/archive-index-data/software/bin/howdesbt cluster --list=${HOW_DIR}/leafnames --bits=12548804 --tree=${HOW_DIR}/union.sbt --nodename=node{number} --keepallnodes
rm ${HOW_DIR}/leafnames

# build index
/project/archive-index-data/software/bin/howdesbt build --HowDe --tree=${HOW_DIR}/union.sbt --outtree=${HOW_DIR}/howde.sbt

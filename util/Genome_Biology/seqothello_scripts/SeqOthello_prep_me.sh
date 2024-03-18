#!/bin/bash

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

source /project/archive-index-data/smehringer/benchmark.variables

OTH_DIR="${WORKDIR}/SeqOthello_bench"

# create data input file
find ${DATA_DIR} -name *.fna* > ${OTH_DIR}/data.txt

# Do jellyfish
#mkdir -p ${OTH_DIR}/kmers
#parallel --line-buffer -j${NUM_THREADS} -a "${OTH_DIR}/data.txt" $(dirname "$0")/call_jellyfish.sh

./example/ConvertToBinary.sh
./example/MakeGroup.sh

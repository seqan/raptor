#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -e
SCRIPT_ROOT=$(dirname $(readlink -f $0))
source $SCRIPT_ROOT/variables.sh

for FASTA in $(seq -f "$INPUT_DIR/$BIN_NUMBER/bins/bin_%0${#BIN_NUMBER}.0f.fasta" 0 1 $((BIN_NUMBER-1))); do
    $HELPER_BINARY --input ${FASTA} \
                   --output $WORKING_DIRECTORY/bins/$(basename ${FASTA} .fasta).fastq
done

awk 'NR%4==2' $INPUT_DIR/$BIN_NUMBER/reads_e$ERRORS\_$READ_LENGTH/all_10.fastq > $READ_FILE



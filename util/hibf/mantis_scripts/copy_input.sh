#!/usr/bin/env bash
# -----------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/master/LICENSE.md
# -----------------------------------------------------------------------------------------------------

set -e
SCRIPT_ROOT=$(dirname $(readlink -f $0))
source $SCRIPT_ROOT/variables.sh

for FASTA in $(seq -f "$INPUT_DIR/$BIN_NUMBER/bins/bin_%0${#BIN_NUMBER}.0f.fasta" 0 1 $((BIN_NUMBER-1))); do
    $HELPER_BINARY --input ${FASTA} \
                   --output $WORKING_DIRECTORY/bins/$(basename ${FASTA} .fasta).fastq
done

awk 'NR%4==2' $INPUT_DIR/$BIN_NUMBER/reads_e$ERRORS\_$READ_LENGTH/all_10.fastq > $READ_FILE



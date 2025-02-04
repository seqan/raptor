#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -Eeuo pipefail
export LC_ALL=C
export LANG=C

## -----------------------------------------------------------------------------
## Input
## -----------------------------------------------------------------------------

# Directory containing fastq files with all 32-mers of each read.
# See create_reads.sh.
READFILEDIR=/dev/shm/seiler/read_files
# The user_bin.ids file. See compare.sh.
USERBINS=/dev/shm/seiler/user_bin.ids
# Directory containing yara binaries.
YARA=/dev/shm/seiler/bin
# Directory containing yara indices. See create_indices.sh.
YARA_INDICIES_DIR=/dev/shm/seiler/yara_indices_for_comparison_raptor_mantis
# Temporary directory to use.
TMPTMP=/dev/shm/smehringer/

# Argument 1 is a line from a file containing FN/FP.
# See validate_on_existing_index.sh
LINE="$1"

## -----------------------------------------------------------------------------
## Output
## -----------------------------------------------------------------------------

# For ${LINE}="readID:binID":
# Via stdout: "${LINE}:<kmer>", where <kmer> is the number of 32-mers of
# read readID that map exactly to binID.

## -----------------------------------------------------------------------------
## Script
## -----------------------------------------------------------------------------

# Split string, e.g. `GCF_000005845.2_ASM584v2_genomic300:2164`,
# into readID and user bin ID.
# IFS=<separator> read var1 var2 <<< <string>
IFS=: read readID binID <<< $LINE

# Extract accession from string, e.g.
# `2164    GCF_000332755.1_ASM33275v1_genomic`. m 1 => Stop after 1 match
accession=$(rg -m 1 "^${binID}" ${USERBINS} | cut -f 2)

# Yara mapping and evaluation.
index=${YARA_INDICIES_DIR}/${accession}/index
alignment=$(mktemp --tmpdir=${TMPTMP} tmp.XXXXXXXXXX.sam)
readkmer_query_file=${READFILEDIR}/${readID}.fastq
${YARA}/yara_mapper \
    ${index} \
    ${readkmer_query_file} \
    -o ${alignment} \
    --threads 1 \
    --error-rate 0 \
    --strata-count 0 \
    --sensitivity full

echo -e "$LINE:$(rg --include-zero -c "NM:i:0" ${alignment})"

rm ${alignment}

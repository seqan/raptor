#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -Eeuo pipefail
export LC_ALL=C
export LANG=C

## -----------------------------------------------------------------------------
## Input
## -----------------------------------------------------------------------------

# Output directory. Will contain files <READID>.fastq with all 32-mers of that read.
OUTPUT_DIR=/dev/shm/seiler/read_files
# The source read file.
READFILE=/dev/shm/seiler/RefSeqCG_arc_bac-queries-1mMio-length250-2errors.fastq.only250.fastq

# Argument 1 is the read ID.
READID="$1"

## -----------------------------------------------------------------------------
## Output
## -----------------------------------------------------------------------------

# ${OUTPUT_DIR}/${READID}.fastq

## -----------------------------------------------------------------------------
## Script
## -----------------------------------------------------------------------------

# Colors
DEFAULT='\033[0m'
RED='\033[0;31m'
GREEN='\033[0;32m'
GRAY='\033[0;90m'

function format_file ()
{
    local file_dirname=$(dirname $1)
    local file_basename=$(basename $1)
    echo -e "${GRAY}${file_dirname}/${DEFAULT}${file_basename}"
}

output_file=${OUTPUT_DIR}/${READID}.fastq

if [[ -s ${output_file} ]]; then
    echo -e "${RED}## Skipping read file generation:${DEFAULT} $(format_file ${output_file})"
else
    echo -e "${GREEN}## Read file generation:${DEFAULT} $(format_file ${output_file})"
    sequence=$(rg -m 1 -A 1 ${READID} ${READFILE} | tail -n 1)
    quality=$(printf %32s | tr \  I) # Produce 32 spaces, then replace space with I.
    for (( i = 0; i+31 < ${#sequence}; i += 1 )); do
        printf '@%s\n%s\n+\n%s\n' "${i}" "${sequence:i:32}" "${quality}" >> ${output_file}
    done
fi

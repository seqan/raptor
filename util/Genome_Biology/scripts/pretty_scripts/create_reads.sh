#!/usr/bin/env bash
set -Eeuo pipefail
export LC_ALL=C
export LANG=C

## -----------------------------------------------------------------------------
## Input
## -----------------------------------------------------------------------------

# Argument 1 is a file with one read ID per line.
# For example, generated via
# cut -f 1 -d: raptor.fns | sort -u > uniq_reads.fns
INPUT_FILE="$1"

## -----------------------------------------------------------------------------
## Output
## -----------------------------------------------------------------------------

# See _create_reads.sh.

## -----------------------------------------------------------------------------
## Script
## -----------------------------------------------------------------------------

if ! [[ -s ${INPUT_FILE} ]]; then
    echo "File ${INPUT_FILE} does not exist or is empty. Exiting."
    exit 1
fi

output_dir=$(grep -o "OUTPUT_DIR=.*" $(dirname "$0")/_create_reads.sh | cut -c 12-)
mkdir -p ${output_dir}
cp "$0" ${output_dir}/$(basename $0).sh
cp "$(dirname "$0")/_create_reads.sh" ${output_dir}/_create_reads.sh

# Call _create_reads.sh in parallel. With 70 threads. Do per line buffering of the stdout/stderr stream.
# Each line in the INPUT_FILE is a argument to call the script with.
parallel --line-buffer -j70 -a "${INPUT_FILE}" $(dirname "$0")/_create_reads.sh

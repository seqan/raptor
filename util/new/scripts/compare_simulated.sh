#!/usr/bin/env bash
set -Eeuo pipefail
LANG=C # For sort

## -----------------------------------------------------------------------------
## Input
## -----------------------------------------------------------------------------

# The original read file.
READFILE=/project/archive-index-data/data/simulated_data/65536/reads_e2_250/all_10.fastq
# Mantis query results file.
MANTIS_RESULTS=/project/archive-index-data/seiler/nature_hibf/runs/mantis/65536/result.txt
# Raptor query results file.
RAPTOR_RESULTS=/dev/shm/seiler/raptor_bench/65536/0.015/24_20_256.out
# Contains compare_mantis_raptor_output and normalise_mantis_output.
RAPTOR_UTIL=/project/archive-index-data/seiler/raptor/build/util2_debug/bin
# Where to write output.
WORKDIR=/project/archive-index-data/seiler/nature_hibf/runs/65K_raptor_mantis

## -----------------------------------------------------------------------------
## Output
## -----------------------------------------------------------------------------

# ${WORKDIR}/${readID}.fastq
# ${WORKDIR}/index_dir/${accession}/index.{lf.{drp,drs,drv,pst}, rid.{concat,limits}, sa.{ind,len,val}, txt.{concat,limits,size}}
# ${WORKDIR}/${readID}_X_${accession}.sam
# Via stdout:
#   * Status messages prefixed with `##`.
#   * Result: "Found \d+ k-mers"

## -----------------------------------------------------------------------------
## Script
## -----------------------------------------------------------------------------

mkdir -p ${WORKDIR}
cp "$0" ${WORKDIR}/script.sh

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

# Extract query IDs from read file. Must be the same order as given to mantis.
query_names=${WORKDIR}/query.names
if [[ -s ${query_names} ]]; then
    echo -e "${RED}## Skipping query name generation: $(format_file ${query_names})"
else
    echo -e "${GREEN}## Query name generation: $(format_file ${query_names})"
    grep '@' $READFILE | tr --delete '@' > ${query_names}
fi

# Extract user bin IDs from raptor results.
user_bin_ids=${WORKDIR}/user_bin.ids
if [[ -s ${user_bin_ids} ]]; then
    echo -e "${RED}## Skipping user bin ID generation: $(format_file ${user_bin_ids})"
else
    echo -e "${GREEN}## User bin ID generation: $(format_file ${user_bin_ids})"
    grep --extended-regexp "^#" ${RAPTOR_RESULTS} \
    | sed 's;/.*/;;g' \
    | sed 's;\.fna.gz;;g' \
    | sed 's;\.fasta;;g' \
    | sed 's;\.minimiser;;g' \
    | tr --delete '#' \
    > ${user_bin_ids}
fi

# Normalise mantis results.
mantis_normalised=${WORKDIR}/mantis_normalised.txt
if [[ -s ${mantis_normalised} ]]; then
    echo -e "${RED}## Skipping normalising mantis output: $(format_file ${mantis_normalised})"
else
    echo -e "${GREEN}## Normalising mantis output: $(format_file ${mantis_normalised})"
    ${RAPTOR_UTIL}/normalise_mantis_output \
        --query_names ${query_names} \
        --user_bin_ids ${user_bin_ids} \
        --mantis_results ${MANTIS_RESULTS} \
        --output_file /dev/stdout \
        --percentage 0.70 \
    2> >(awk '$0="    "$0' >&2) \
    | sort > ${mantis_normalised}
fi

# Normalise raptor results.
raptor_normalised=${WORKDIR}/raptor_normalised.txt
if [[ -s ${raptor_normalised} ]]; then
    echo -e "${RED}## Skipping normalising raptor output: $(format_file ${raptor_normalised})"
else
    echo -e "${GREEN}## Normalising raptor output: $(format_file ${raptor_normalised})"
    grep -v "^#" ${RAPTOR_RESULTS} | sort > ${raptor_normalised}
fi

# Compare results.
echo -e "${GREEN}## Comparing $(format_file ${mantis_normalised})${GREEN} and $(format_file ${raptor_normalised})"
${RAPTOR_UTIL}/compare_mantis_raptor_output \
    --mantis_results ${mantis_normalised} \
    --raptor_results ${raptor_normalised} \
    --user_bin_ids ${user_bin_ids} \
    --output_directory ${WORKDIR} \
    | awk '$0="    "$0'

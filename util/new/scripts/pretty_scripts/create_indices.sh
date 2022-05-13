#!/usr/bin/env bash
set -Eeuo pipefail
export LC_ALL=C
export LANG=C

## -----------------------------------------------------------------------------
## Input
## -----------------------------------------------------------------------------

# Directory containing all RefSeq data (.../v1/files).
REFSEQ=/project/archive-index-data/data/RefSeqCG_arc_bac/v1/files
# The user_bin.ids file. See compare.sh.
USERBINS=/project/archive-index-data/seiler/kmer_RefSeq/compare/user_bin.ids
# Directory containing yara binaries.
YARA=/project/archive-index-data/seiler/kmer_RefSeq/yara/build/bin
# Where to write output indices to.
OUTPUT_DIR=/project/archive-index-data/smehringer/yara_indices_for_comparison_raptor_mantis/

## -----------------------------------------------------------------------------
## Output
## -----------------------------------------------------------------------------

# A yara index for each accession in ${USERBINS}:
# ${OUTPUT_DIR}/${accession}/index.{lf.{drp,drs,drv,pst}, rid.{concat,limits}, sa.{ind,len,val}, txt.{concat,limits,size}}

## -----------------------------------------------------------------------------
## Script
## -----------------------------------------------------------------------------

mkdir -p ${OUTPUT_DIR}
cp "$0" ${OUTPUT_DIR}/$(basename $0).sh

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

for accession in $(cut -f2 ${USERBINS})
do
    index=${OUTPUT_DIR}/${accession}/index
    if [[ $(ls -1 ${index}* 2>/dev/null | wc -l ) -eq 12 ]]; then
        echo -e "${RED}## Skipping index generation:${DEFAULT} $(format_file ${index})* already exist."
    else
        echo -e "${GREEN}## Index generation:${DEFAULT} $(format_file ${index})*"
        mkdir -p $(dirname ${index})
        ${YARA}/yara_indexer \
            ${REFSEQ}/${accession}.fna.gz \
            -o ${index}
    fi
done

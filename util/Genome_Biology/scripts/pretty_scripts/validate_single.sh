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

# Contains user_bins.ids and raptor.fns.
COMPAREDIR=/project/archive-index-data/seiler/kmer_RefSeq/compare
# The original read file.
READFILE=/project/archive-index-data/data/RefSeqCG_arc_bac/RefSeqCG_arc_bac-queries-1mMio-length250-2errors.fastq
# Contains RefSeq sequence files, i.e. <accession>.fna.gz.
REFSEQ=/project/archive-index-data/data/RefSeqCG_arc_bac/v1/files
# Contains yara_indexer and yara_mapper.
YARA=/project/archive-index-data/seiler/kmer_RefSeq/yara/build/bin
# Where to write output.
WORKDIR=/project/archive-index-data/seiler/kmer_RefSeq/yara/validation

## -----------------------------------------------------------------------------
## Output
## -----------------------------------------------------------------------------

# ${WORKDIR}/script.sh
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
cp "$0" ${WORKDIR}/$(basename $0).sh

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

# Split string, e.g. `GCF_000005845.2_ASM584v2_genomic300:2164`, into readID and user bin ID.
# IFS=<separator> read var1 var2 <<< <string>
# IFS=: read readID binID <<< $(head -n 1 ${COMPAREDIR}/compare/raptor.fns)
readID=GCF_000757015.2_ASM75701v2_genomic370
binID=4068

# Extract accession from string, e.g. `2164    GCF_000332755.1_ASM33275v1_genomic`. m 1 => Stop after 1 match
accession=$(grep -m 1 "^${binID}" ${COMPAREDIR}/user_bin.ids | cut -f 2)

# Generate kmers to search.
readkmer=${WORKDIR}/${readID}.fastq
if [[ -s ${readkmer} ]]; then
    echo -e "${RED}## Skipping k-mer generation:${DEFAULT} $(format_file ${readkmer})"
else
    echo -e "${GREEN}## k-mer generation:${DEFAULT} $(format_file ${readkmer})"
    # Get read sequence. m 1 => Stop after one match, A 1 => Also print one line after match.
    # For example:
    #  @GCF_000005845.2_ASM584v2_genomic300
    #  GTTGCAAACTGGTGC[...]
    # tail then takes last line only
    sequence=$(grep -m 1 -A 1 ${readID} ${READFILE} | tail -n 1)
    quality=$(printf %32s | tr \  I) # Produce 32 spaces, then replace space with I.
    for (( i = 0; i+31 < ${#sequence}; i += 1 )); do
        printf '@%s\n%s\n+\n%s\n' "${i}" "${sequence:i:32}" "${quality}" >> ${readkmer}
    done
fi

# Yara index.
index=${WORKDIR}/index_dir/${accession}/index
if [[ $(ls -1 ${index}* 2>/dev/null | wc -l ) -eq 12 ]]; then
    echo -e "${RED}## Skipping index generation:${DEFAULT} $(format_file ${index})* already exist."
else
    echo -e "${GREEN}## Index generation:${DEFAULT} $(format_file ${index})*"
    mkdir -p $(dirname ${index})
    ${YARA}/yara_indexer \
        ${REFSEQ}/${accession}.fna.gz \
        -o ${index}
fi

# Yara mapping.
alignment=${WORKDIR}/${readID}_X_${accession}.sam
if [[ -s ${alignment} ]]; then
    echo -e "${RED}## Skipping alignment generation:${DEFAULT} $(format_file ${alignment})"
else
    echo -e "${GREEN}## Alignment generation:${DEFAULT} $(format_file ${alignment})"
    ${YARA}/yara_mapper \
        ${index} \
        ${readkmer} \
        -o ${alignment} \
        --error-rate 0 \
        --strata-count 0 \
        --secondary-matches record \
        --align-secondary \
        --sensitivity full
fi

echo -e "Found \e[1m$(grep -c "NM:i:0" ${alignment})\e[0m k-mers"

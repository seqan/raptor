#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -eEuo pipefail

DATA_DIR="/project/archive-index-data/data/RefSeqCG_arc_bac/v1/files/"
DATA_FILENAMES="/project/archive-index-data/data/RefSeqCG_arc_bac/filenames.txt"

QUERY_1K="/project/archive-index-data/data/RefSeqCG_arc_bac/RefSeqCG_arc_bac-queries-1T-length250-2errors.fastq.only250.fastq"
QUERY_1M="/project/archive-index-data/data/RefSeqCG_arc_bac/RefSeqCG_arc_bac-queries-1M-length250-2errors.fastq.only250.fastq"
QUERY_5M="/project/archive-index-data/data/RefSeqCG_arc_bac/RefSeqCG_arc_bac-queries-5M-length250-2errors.fastq.only250.fastq"
QUERY_10M="/project/archive-index-data/data/RefSeqCG_arc_bac/RefSeqCG_arc_bac-queries-1mMio-length250-2errors.fastq.only250.fastq"

WORK_DIR="/project/archive-index-data/seiler/refseq_metagraph"
LOG_DIR="${WORK_DIR}/logs"
TIME_DIR="${WORK_DIR}/times"
KMC_DIR="${WORK_DIR}/kmc_files"
INDEX_DIR="${WORK_DIR}/index"
ANNO_DIR="${INDEX_DIR}/annotations" # separate annotations

mkdir -p ${LOG_DIR}
mkdir -p ${TIME_DIR}
mkdir -p ${KMC_DIR}
mkdir -p ${INDEX_DIR}
mkdir -p ${ANNO_DIR}

KMC_FILENAMES="${KMC_DIR}/file_list.txt"

KMER_SIZE=32
NUM_THREADS=32
QUERY_THRESHOLD=0.7

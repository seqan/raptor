#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -Eeuo pipefail

echo -n "[$(date +"%Y-%m-%d %T")] Running KMC..."
/project/archive-index-data/seiler/refseq_metagraph/0_kmc.sh
echo "Done."

echo -n "[$(date +"%Y-%m-%d %T")] Running first build..."
/project/archive-index-data/seiler/refseq_metagraph/1_build.first.sh
echo "Done."

echo -n "[$(date +"%Y-%m-%d %T")] Running transform..."
/project/archive-index-data/seiler/refseq_metagraph/2_transform.sh
echo "Done."

echo -n "[$(date +"%Y-%m-%d %T")] Running second build..."
/project/archive-index-data/seiler/refseq_metagraph/3_build.second.sh
echo "Done."

echo -n "[$(date +"%Y-%m-%d %T")] Running annotate..."
/project/archive-index-data/seiler/refseq_metagraph/4_annotate.sh
echo "Done."

echo -n "[$(date +"%Y-%m-%d %T")] Running transform_anno..."
/project/archive-index-data/seiler/refseq_metagraph/5_transform_anno.sh
echo "Done."

echo -n "[$(date +"%Y-%m-%d %T")] Running query..."
/project/archive-index-data/seiler/refseq_metagraph/6_query.sh
echo "Done."

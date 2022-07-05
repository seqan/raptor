#!/usr/bin/env bash
# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------

set -e
SCRIPT_ROOT=$(dirname $(readlink -f $0))
source $SCRIPT_ROOT/variables.sh

trap 'echo; echo "## [$(date +"%Y-%m-%d %T")] ERROR ##"' ERR

echo "## [$(date +"%Y-%m-%d %T")] Start ##"

echo -n "[$(date +"%Y-%m-%d %T")] Copying input..."
    /usr/bin/time -o $copy_time -v \
        $SCRIPT_ROOT/copy_input.sh \
            &>$copy_log
echo "Done."

echo -n "[$(date +"%Y-%m-%d %T")] Running squeakr..."
    /usr/bin/time -o $squeakr_time -v \
        $SCRIPT_ROOT/squeakr.sh \
            &>$squeakr_log
echo "Done."

echo -n "[$(date +"%Y-%m-%d %T")] Running mantis build..."
    /usr/bin/time -o $mantis_build_time -v \
        $SCRIPT_ROOT/mantis_build.sh \
            &>$mantis_build_log
echo "Done."

echo -n "[$(date +"%Y-%m-%d %T")] Running mantis mst..."
    /usr/bin/time -o $mantis_mst_time -v \
        $SCRIPT_ROOT/mantis_mst.sh \
            &>$mantis_mst_log
echo "Done."

echo -n "[$(date +"%Y-%m-%d %T")] Running mantis query..."
    /usr/bin/time -o $mantis_query_time -v \
        $SCRIPT_ROOT/mantis_query.sh \
            &>$mantis_query_log
echo "Done."

echo "## [$(date +"%Y-%m-%d %T")] End ##"

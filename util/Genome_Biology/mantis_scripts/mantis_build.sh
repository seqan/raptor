#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

set -e
SCRIPT_ROOT=$(dirname $(readlink -f $0))
source $SCRIPT_ROOT/variables.sh

seq -f "$SQUEAKR_DIRECTORY/bin_%0${#BIN_NUMBER}.0f.squeakr" 0 1 $((BIN_NUMBER-1)) > $WORKING_DIRECTORY/bin_list.tmp
ulimit -Sn $(ulimit -Hn)
$MANTIS_BINARY build \
    -s 34 \
    -i $WORKING_DIRECTORY/bin_list.tmp \
    -o $MANTIS_INDEX
rm $WORKING_DIRECTORY/bin_list.tmp

#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2024 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2024 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

GCC="/usr/bin/g++-11"
DO_TIME=0

if [[ DO_TIME -eq 0 ]]; then
    exec "$GCC" "$@"
else
    FILE=$(mktemp ram_usage.XXXXXXXX)
    exec /usr/bin/time -v "$GCC" "$@" 2> $FILE
fi

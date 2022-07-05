#!/usr/bin/env bash
# --------------------------------------------------------------------------------------------------
# Copyright (c) 2006-2022, Knut Reinert & Freie Universität Berlin
# Copyright (c) 2016-2022, Knut Reinert & MPI für molekulare Genetik
# This file may be used, modified and/or redistributed under the terms of the 3-clause BSD-License
# shipped with this file and also available at: https://github.com/seqan/raptor/blob/main/LICENSE.md
# --------------------------------------------------------------------------------------------------
#
# Usage: collect_compile_stats.sh
# Will compile the project in the current working directory and create a summary of compile time and memory usage per
# compiled binary.
#
# Output will be written to current working directory.

set -e

# Get the full path to this bash script.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"

# We change a variable in gcc.sh and g++.sh, which we will need to reset.
reset_scripts() {
    sed -i "s/DO_TIME=1/DO_TIME=0/" $SCRIPT_DIR/gcc.sh
    sed -i "s/DO_TIME=1/DO_TIME=0/" $SCRIPT_DIR/g++.sh
}
# reset_scripts() will be executed whenever this script exits, e.g. on completion or on error.
trap reset_scripts EXIT

# Run in a subshell `(...)` to `set -x` just for one command.
(set -ex; cmake $SCRIPT_DIR/../.. -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_COMPILER=$SCRIPT_DIR/g++.sh -DCMAKE_C_COMPILER=$SCRIPT_DIR/gcc.sh -DRAPTOR_USE_CCACHE=OFF -DRAPTOR_NATIVE_BUILD=ON)

# We need `DO_TIME=0` for the CMake configuration to work, but `DO_TIME=1` to actually measure time/RAM consumption.
sed -i "s/DO_TIME=0/DO_TIME=1/" $SCRIPT_DIR/gcc.sh
sed -i "s/DO_TIME=0/DO_TIME=1/" $SCRIPT_DIR/g++.sh

# Run in a subshell `(...)` to `set -x` just for one command.
(set -ex; make -k -j6 cli_test api_test)

# Concatenate results and run python script for fomatted output.
find . -name "ram_usage.*" -exec cat {} + > complete.txt
$SCRIPT_DIR/ram_usage.py complete.txt stats.csv

echo "Results can be found in $(pwd)/stats.csv"

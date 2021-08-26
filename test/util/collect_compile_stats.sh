#!/usr/bin/env bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"

reset_scripts() {
    sed -i "s/DO_TIME=1/DO_TIME=0/" $SCRIPT_DIR/gcc.sh
    sed -i "s/DO_TIME=1/DO_TIME=0/" $SCRIPT_DIR/g++.sh
}
trap reset_scripts EXIT

set -ex

cmake $SCRIPT_DIR/../.. -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_COMPILER=$SCRIPT_DIR/g++.sh -DCMAKE_C_COMPILER=$SCRIPT_DIR/gcc.sh -DUSE_CCACHE=OFF

sed -i "s/DO_TIME=0/DO_TIME=1/" $SCRIPT_DIR/gcc.sh
sed -i "s/DO_TIME=0/DO_TIME=1/" $SCRIPT_DIR/g++.sh

make -k -j4 cli_test api_test

find . -name "ram_usage.*" -exec cat {} + > complete.txt
$SCRIPT_DIR/parse.py complete.txt stats.csv

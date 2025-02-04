#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2006-2025 Knut Reinert & Freie Universität Berlin
# SPDX-FileCopyrightText: 2016-2025 Knut Reinert & MPI für molekulare Genetik
# SPDX-License-Identifier: BSD-3-Clause

cutoff=50
cutoffs=(1 3 10 20)
size=0
sizes=(314572800 524288000 1073741824 3221225472)
for file in "$@"
do
    size=$(du -b $file | awk '{ print $1 }')
    echo $size
    basename "$file"
    name="$(basename -- $file)"
    for i in 0 1 2 3
    do
        if (($size < ${sizes[$i]}));
        then
        echo ${sizes[$i]}
        cutoff=${cutoffs[$i]}
        break
        fi
    done
    echo $name
    echo $cutoff
    ./squeakr count -e -k 20 -n -c $cutoff -t 1 -o ${name}.squeakr $file
    cutoff=50
done

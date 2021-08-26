#!/usr/bin/env bash

GCC="/usr/bin/g++-11"
DO_TIME=0

if [[ DO_TIME -eq 0 ]]; then
    exec "$GCC" "$@"
else
    FILE=$(mktemp ram_usage.XXXXXXXX)
    exec /usr/bin/time -v "$GCC" "$@" 2> $FILE
fi

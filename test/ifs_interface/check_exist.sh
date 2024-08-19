#!/usr/bin/env bash

set -e

FILE=./lib/libfesom.a
if [ -f "$FILE" ]; then
    echo "$FILE compiled and linked."
    exit 0
else
    echo "$FILE does not exist."
    exit 1
fi

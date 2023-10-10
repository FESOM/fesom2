#!/usr/bin/env bash

set -e

LIB_NAMES=("libfesom.a" "libfesom.so" "libfesom.dylib")

libfesom_exists=false
for LIB_NAME in ${LIB_NAMES[@]}; do
  FILE=./lib/${LIB_NAME}
  if [ -f "$FILE" ]; then
    echo "$FILE compiled and linked."
    libfesom_exists=true
  fi
done

if ${libfesom_exists}; then
    exit 0
else
    echo "libfesom not found"
    exit 1
fi


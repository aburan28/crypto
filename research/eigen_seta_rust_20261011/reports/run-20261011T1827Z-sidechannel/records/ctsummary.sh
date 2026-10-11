#!/bin/sh
# usage: ctsummary.sh DIR  -> table of total / libc allocator+mem* / rest per run
cd "$1"
printf "%-16s %15s %12s %15s\n" run total alloc_mem rest
for f in cg.*.*.*; do
  case $f in *.fn|*.sg) continue;; esac
  callgrind_annotate --auto=no --threshold=100 $f 2>/dev/null | grep -E "^ *[0-9,]+ \(" | sed -E 's/^ *([0-9,]+) \([^)]*\) +/\1\t/' | awk -F'\t' -v name=$f '
    {gsub(",","",$1)}
    $2 ~ /PROGRAM TOTALS/ {tot=$1; next}
    $2 ~ /malloc\/|arena\.c|memset|memmove|memcpy|string_fortified|\?\?\?:0x/ {al+=$1}
    END {printf "%-16s %15d %12d %15d\n", name, tot, al, tot-al}'
done

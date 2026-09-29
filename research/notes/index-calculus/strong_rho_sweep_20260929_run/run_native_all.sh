#!/bin/bash
cd "$(dirname "$0")"
while read -r N L K CORPUS; do
  echo "== native n=$N L=$L K=$K ($(date +%T))"
  ./sweep_cell.py native "$N" "$L" "$K" "$CORPUS" 2>&1 | cut -c1-600
done < cells.txt
echo NATIVE_ALL_DONE

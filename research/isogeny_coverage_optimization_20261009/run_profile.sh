#!/usr/bin/env bash
# Thin orchestration of the frozen native binaries and Valgrind.
set -euo pipefail
mkdir -p /evidence
for case_name in screen-p224 verify-p224 search-p192; do
  for round in 1 2 3 4 5; do
    for build in baseline candidate; do
      run_dir="/evidence/$case_name-$round-$build"
      mkdir -p "$run_dir"
      bin="/build/$build/release/isogeny"
      case "$case_name" in
        screen-p224) args=(screen --curve p224 --from 1010 --to 5000 --max-order 8) ;;
        verify-p224) args=(verify --input /fixtures/p224-1471) ;;
        search-p192) args=(search --curve p192 --ell 73 --out "$run_dir/results" --timeout 1800) ;;
      esac
      printf '%s %s %s\n' "$case_name" "$round" "$build"
      timeout 1800 valgrind --tool=callgrind --trace-children=yes \
        --callgrind-out-file="$run_dir/callgrind.%p" \
        "$bin" "${args[@]}" > "$run_dir/stdout.json" 2> "$run_dir/stderr.log"
      awk '/^summary:/ {total += $2} END {printf "%.0f\n", total}' "$run_dir"/callgrind.* > "$run_dir/instructions.txt"
      cat "$run_dir/instructions.txt"
    done
  done
done

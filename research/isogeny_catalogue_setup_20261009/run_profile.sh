#!/usr/bin/env bash
# Thin orchestration only: algorithms, verification and analysis are native Rust.
set -euo pipefail
mkdir -p /evidence
rustc -vV > /evidence/rustc.txt
valgrind --version > /evidence/valgrind.txt
uname -a > /evidence/kernel.txt
for build in baseline candidate; do
  /build/$build/release/isogeny --version > "/evidence/$build-version.json"
  sha256sum "/build/$build/release/isogeny" > "/evidence/$build-binary.sha256"
  /build/$build/release/isogeny curves > "/evidence/$build-catalogue.json"
done
for case_name in screen-p224 verify-p224 search-p192; do
  for round in 1 2 3 4 5; do
    for build in baseline candidate; do
      run_dir="/evidence/$case_name-$round-$build"
      mkdir -p "$run_dir"
      case "$case_name" in
        screen-p224) args=(screen --curve p224 --from 1010 --to 5000 --max-order 8) ;;
        verify-p224) args=(verify --input /fixtures/p224-1471) ;;
        search-p192) args=(search --curve p192 --ell 73 --out "$run_dir/results" --timeout 1800) ;;
      esac
      printf '%s %s %s\n' "$case_name" "$round" "$build"
      set +e
      timeout 1800 valgrind --tool=callgrind --trace-children=yes \
        --callgrind-out-file="$run_dir/callgrind.%p" \
        "/build/$build/release/isogeny" "${args[@]}" \
        > "$run_dir/stdout.json" 2> "$run_dir/stderr.log"
      exit_status=$?
      set -e
      printf '%s\n' "$exit_status" > "$run_dir/exit-status.txt"
      if test "$exit_status" -ne 0; then
        printf 'FAILED: %s (exit %s)\n' "$run_dir" "$exit_status" >&2
        exit "$exit_status"
      fi
      awk '/^summary:/ {total += $2} END {printf "%.0f\n", total}' \
        "$run_dir"/callgrind.* > "$run_dir/instructions.txt"
      cat "$run_dir/instructions.txt"
    done
  done
done

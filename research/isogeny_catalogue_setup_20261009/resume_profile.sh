#!/usr/bin/env bash
# Resume thin native orchestration without rewriting completed or failed records.
set -euo pipefail
mkdir -p /evidence/environment-after-restart
rustc -vV > /evidence/environment-after-restart/rustc.txt
valgrind --version > /evidence/environment-after-restart/valgrind.txt
uname -a > /evidence/environment-after-restart/kernel.txt
cat /proc/cpuinfo > /evidence/environment-after-restart/cpuinfo.txt
cmp /evidence/rustc.txt /evidence/environment-after-restart/rustc.txt
cmp /evidence/valgrind.txt /evidence/environment-after-restart/valgrind.txt
for build in baseline candidate; do
  /build/$build/release/isogeny --version > "/evidence/environment-after-restart/$build-version.json"
  sha256sum "/build/$build/release/isogeny" > "/evidence/environment-after-restart/$build-binary.sha256"
  cmp "/evidence/$build-version.json" "/evidence/environment-after-restart/$build-version.json"
  cmp "/evidence/$build-binary.sha256" "/evidence/environment-after-restart/$build-binary.sha256"
done
for case_name in screen-p224 verify-p224 search-p192; do
  for round in 1 2 3 4 5; do
    for build in baseline candidate; do
      run_dir="/evidence/$case_name-$round-$build"
      if test -f "$run_dir/exit-status.txt" && test -s "$run_dir/instructions.txt"; then
        test "$(cat "$run_dir/exit-status.txt")" = 0
        printf 'RETAINED %s %s %s\n' "$case_name" "$round" "$build"
        continue
      fi
      test ! -e "$run_dir" || { printf 'Unsealed run must be preserved separately: %s\n' "$run_dir" >&2; exit 2; }
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
        "/build/$build/release/isogeny" "${args[@]}" > "$run_dir/stdout.json" 2> "$run_dir/stderr.log"
      exit_status=$?
      set -e
      printf '%s\n' "$exit_status" > "$run_dir/exit-status.txt"
      test "$exit_status" -eq 0 || exit "$exit_status"
      awk '/^summary:/ {total += $2} END {printf "%.0f\n", total}' "$run_dir"/callgrind.* > "$run_dir/instructions.txt"
      cat "$run_dir/instructions.txt"
    done
  done
done

#!/usr/bin/env bash
# Thin Linux orchestration: the paired benchmark and analysis are native Rust.
set -u

if [ "$#" -ne 4 ]; then
  echo 'usage: run_linux.sh THREADS RESERVE_CPUS PIN_CPUS OUTPUT_DIR' >&2
  exit 2
fi
threads="$1"
reserve="$2"
pin="$3"
root="$4"
lock="$RUNNER_TEMP/f5-support-split-bench.lock"
overall=0
for seed in frozen holdout_a holdout_b holdout_c; do
  selected=0
  for attempt in 1 2 3; do
    directory="$root/$seed/attempt$attempt"
    mkdir -p "$directory"
    timeout 180s sudo -n env \
      "PATH=$PATH" "HOME=$HOME" \
      "RUSTUP_HOME=${RUSTUP_HOME:-$HOME/.rustup}" \
      "CARGO_HOME=${CARGO_HOME:-$HOME/.cargo}" \
      "GITHUB_HEAD_REF=${GITHUB_HEAD_REF:-}" \
      python3 tools/isolated_bench.py reserve \
      --lock "$lock" --cpus "$reserve" \
      --out "$directory/isolation.jsonl" --settle 10 --period 2 \
      --label "f5-support-split-t$threads-$seed-attempt$attempt" \
      -- taskset -c "$pin" target/release/examples/f5_support_split_pair \
      target/release/examples/f4_f2_bench "$directory/paired.json" \
      "$seed" "$threads" >"$directory/stdout.txt" 2>"$directory/stderr.txt"
    code=$?
    printf '{"exit_code":%s}\n' "$code" >"$directory/exit.json"
    if [ "$code" -eq 0 ] \
      && tail -n 1 "$directory/isolation.jsonl" | jq -e '.exit_status == 0 and .contended_samples == 0 and ((.left_on_reserved.user_threads | length) == 0)' >/dev/null \
      && jq -e '.status == "complete" and (.runs | length) == 22 and all(.runs[]; .status == "ok")' "$directory/paired.json" >/dev/null; then
      selected="$attempt"
      break
    fi
  done
  printf '{"seed":"%s","selected_attempt":%s}\n' "$seed" "$selected" >"$root/$seed/selection.json"
  if [ "$selected" -eq 0 ]; then overall=1; fi
done
exit "$overall"

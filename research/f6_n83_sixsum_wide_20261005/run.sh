#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 || $# -gt 3 ]]; then
  echo 'usage: run.sh LABEL planted | LABEL ordinary OFFSET' >&2
  exit 2
fi
label=$1
mode=$2
offset=${3-}
case $label in
  *[!a-zA-Z0-9_-]*|'') echo 'invalid label' >&2; exit 2;;
esac
root=$(cd "$(dirname "$0")/../.." && pwd)
dir="$root/research/f6_n83_sixsum_wide_20261005"
binary=/Volumes/SSD990/crypto/worktrees/f6-n83-geometric-20261004/target/release/examples/f6_n83_sixsum_wide_probe
out="$dir/$label.jsonl"
err="$dir/$label.stderr.txt"
status="$dir/$label.status"
limit_kib=$((7 * 1024 * 1024))
seconds_limit=120
max_kib=0
reason=completed

args=("$mode")
if [[ -n $offset ]]; then args+=("$offset"); fi
"$binary" "${args[@]}" >"$out" 2>"$err" &
pid=$!
start=$SECONDS
while kill -0 "$pid" 2>/dev/null; do
  # The process can exit between kill -0 and ps. Treat that race as an
  # empty sample; wait below still records its real exit status.
  rss=$(ps -o rss= -p "$pid" 2>/dev/null | tr -d '[:space:]' || true)
  if [[ $rss =~ ^[0-9]+$ ]] && (( rss > max_kib )); then
    max_kib=$rss
  fi
  if (( max_kib >= limit_kib )); then
    reason=memory_guard
    kill -TERM "$pid" 2>/dev/null || true
    break
  fi
  if (( SECONDS - start >= seconds_limit )); then
    reason=time_guard
    kill -TERM "$pid" 2>/dev/null || true
    break
  fi
  sleep 0.05
done
set +e
wait "$pid"
exit_code=$?
set -e
printf 'reason=%s\nexit_code=%s\nmax_sampled_rss_kib=%s\nelapsed_s=%s\n' \
  "$reason" "$exit_code" "$max_kib" "$((SECONDS-start))" >"$status"
cat "$status"
exit "$exit_code"

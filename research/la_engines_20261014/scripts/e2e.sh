#!/usr/bin/env bash
# End to end (PROTOCOL.md question 3): real `ic logs` precomputes, reading
# elapsed_seconds and linear_algebra.seconds from the report (relation
# collection + certification = their difference, labelled derived).
# usage: e2e.sh  (writes raw/e2e.jsonl)
set -uo pipefail
here=$(cd "$(dirname "$0")/.." && pwd)
IC=${IC_BIN:-/home/user/crypto/target/release/ic}
out=$here/raw/e2e.jsonl
work=${TMPDIR:-/tmp/cb}/e2e
mkdir -p "$work"
TIMEOUT=1800

# Factor bases chosen by `ic search` (target-blind), saved once.
spec() { # degree summands curve_a
  local f="$work/fb$1-m$2-a$3.json"
  [[ -f $f ]] || "$IC" search --degree "$1" --curve-a "$3" --summands "$2" --spec-out "$f" > "$work/search$1-m$2.log" 2>&1
  echo "$f"
}

run() { # label degree summands curve_a la sparse_solver seed [fb]
  local label=$1 n=$2 m=$3 a=$4 la=$5 ss=$6 seed=$7 fb=${8:-}
  local db="$work/db-$label-$la-$ss-$seed-$RANDOM.json"
  local args=(logs --degree "$n" --curve-a "$a" --summands "$m" --seed "$seed"
    --linear-algebra "$la" --database "$db" --json)
  [[ $la == sparse ]] && args+=(--sparse-solver "$ss")
  [[ -n $fb ]] && args+=(--factor-base "$fb")
  local start end rc rep
  start=$(date +%s.%N)
  rep=$(timeout $TIMEOUT "$IC" "${args[@]}" 2> "$work/err-$label-$la-$ss-$seed.log"); rc=$?
  end=$(date +%s.%N)
  local el; el=$(awk -v a="$start" -v b="$end" 'BEGIN{print b-a}')
  if [[ $rc -ne 0 ]]; then
    jq -nc --arg label "$label" --arg la "$la" --arg ss "$ss" --argjson seed "$seed" --argjson rc "$rc" --argjson el "$el" \
      '{label:$label, linear_algebra:$la, sparse_solver:$ss, seed:$seed, failed:true, exit:$rc, process_elapsed_s:$el}' >> "$out"
  else
    echo "$rep" | jq -c --arg label "$label" --arg la "$la" --arg ss "$ss" --argjson seed "$seed" --argjson el "$el" \
      --argjson n "$n" --argjson m "$m" \
      '{label:$label, degree:$n, summands:$m, linear_algebra:$la, sparse_solver:$ss, seed:$seed,
          columns:.counts.columns, relations:.counts.relations, trials:.counts.trials,
          verified:.verified, status:.status, elapsed_seconds:.elapsed_seconds,
          linear_algebra_seconds:.linear_algebra.seconds, la_attempts:.linear_algebra.attempts,
          core_dimension:(.linear_algebra.sparse.core_dimension // null),
          other_seconds_derived:(.elapsed_seconds - .linear_algebra.seconds),
          cpu_seconds:.resources.cpu_seconds, peak_rss_bytes:.resources.peak_rss_bytes,
          process_elapsed_s:$el}' | head -1 >> "$out"
  fi
  echo "e2e $label $la/$ss seed=$seed rc=$rc ${el}s" >&2
}

configs=()
for cfg in "$@"; do configs+=("$cfg"); done
for cfg in "${configs[@]}"; do
  IFS=: read -r label n m a usefb <<< "$cfg"
  fb=""; [[ $usefb == search ]] && fb=$(spec "$n" "$m" "$a")
  for seed in 1 2 3; do
    for mode in sparse:block-wiedemann sparse:block-lanczos dense:-; do
      run "$label" "$n" "$m" "$a" "${mode%%:*}" "${mode##*:}" "$seed" "$fb"
    done
  done
done
echo done >&2

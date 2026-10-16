#!/usr/bin/env bash
# Amendment 1 (PROTOCOL.md §7): the registered systems, dumped by the harness
# from the same seed and target streams, solved to a grevlex Gröbner basis by
# two external engines. Thin orchestration only: one engine process per
# system, one at a time, under a 6 GB address-space cap and a wall cap.
#
#   BIN=... MSOLVE=... SINGULAR=... run_ext.sh <out-dir>
set -u
OUT=${1:?out dir}
BIN=${BIN:?harness binary}
MSOLVE=${MSOLVE:?msolve binary}
SINGULAR=${SINGULAR:?Singular binary}
SEED=20261016
P=65521
CAP=6000000000
WALL=${WALL:-3600}
mkdir -p "$OUT/systems"
{
  echo "commit $(git rev-parse HEAD)"
  echo "harness_sha256 $(sha256sum "$BIN" | cut -d' ' -f1)"
  echo "msolve $("$MSOLVE" -V 2>&1 | head -1) sha256 $(sha256sum "$MSOLVE" | cut -d' ' -f1)"
  echo "singular $("$SINGULAR" -v --version 2>&1 | head -1) sha256 $(sha256sum "$SINGULAR" | cut -d' ' -f1)"
  echo "host $(uname -n) $(nproc) cores $(free -g | awk '/Mem:/{print $2}') GB"
  echo "start $(date -u +%Y-%m-%dT%H:%M:%SZ)"
} > "$OUT/manifest.txt"

for cell in "2 2 8" "3 2 8" "3 3 8" "4 3 8" "5 4 4" "4 4 4"; do
  set -- $cell
  for kind in planted random; do
    "$BIN" --n "$1" --m "$2" --p $P --seed $SEED --targets "$3" --kind "$kind" \
      --arms sem,rr --dump-dir "$OUT/systems" --dump-only > /dev/null
  done
done
(cd "$OUT/systems" && sha256sum *.ms *.sing > ../systems.sha256)

row() { # tag engine rc wall_s rss_kb extra-json
  local tag=$1
  local n=${tag#n}; n=${n%%m*}
  local m=${tag#*m}; m=${m%%_*}
  local rest=${tag#*_}; local kind=${rest%%_*}; rest=${rest#*_}
  local t=${rest%%_*}; local arm=${rest#*_}
  jq -cn --argjson n "$n" --argjson m "$m" --arg kind "$kind" --argjson t "$t" --arg arm "$arm" \
    --arg engine "$2" --argjson rc "$3" --argjson wall "$4" --argjson rss "$5" --argjson x "$6" \
    '{row:"ext",n:$n,m:$m,kind:$kind,target:$t,arm:$arm,engine:$engine,exit:$rc,wall_s:$wall,rss_kb:$rss} + $x' \
    >> "$OUT/rows.jsonl"
}

for f in $(ls "$OUT"/systems/*.ms | sort -V); do
  tag=$(basename "$f" .ms)
  # msolve: F4 (Faugère's own design lineage), reduced grevlex basis only.
  log=$(mktemp)
  /usr/bin/time -f "TIME %e %M" prlimit --as=$CAP timeout "$WALL" \
    "$MSOLVE" -g 2 -t 1 -v 2 -f "$f" -o /dev/null > "$log" 2>&1
  rc=$?
  read -r wall rss < <(awk '/^TIME/{print $2, $3}' "$log")
  x=$(awk '
    /msolve overall time/ {t=$4}
    /#zero reductions/ {z=$3}
    /#rows reduced/ {r=$3}
    /max. matrix data/ {mr=$4; mc=$6}
    /length of basis/ {b=$5}
    END {printf "{\"engine_time_s\":%s,\"zero_reductions\":%s,\"rows_reduced\":%s,\"max_rows\":%s,\"max_cols\":%s}",
         (t==""?"null":t), (z==""?"null":z), (r==""?"null":r), (mr==""?"null":mr), (mc==""?"null":mc)}' "$log")
  row "$tag" msolve "$rc" "${wall:-null}" "${rss:-null}" "$x"
  rm -f "$log"
  # Singular sba: signature-based Gröbner basis (the F5 family).
  log=$(mktemp)
  /usr/bin/time -f "TIME %e %M" prlimit --as=$CAP timeout "$WALL" \
    "$SINGULAR" -q "${f%.ms}.sing" > "$log" 2>&1
  rc=$?
  read -r wall rss < <(awk '/^TIME/{print $2, $3}' "$log")
  x=$(awk '
    /^sba_ms/ {s=$2} /^vdim/ {v=$2} /^gb_size/ {g=$2}
    END {printf "{\"engine_time_s\":%s,\"vdim\":%s,\"gb_size\":%s}",
         (s==""?"null":s/1000), (v==""?"null":v), (g==""?"null":g)}' "$log")
  row "$tag" sba "$rc" "${wall:-null}" "${rss:-null}" "$x"
  rm -f "$log"
done
echo "end $(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$OUT/manifest.txt"

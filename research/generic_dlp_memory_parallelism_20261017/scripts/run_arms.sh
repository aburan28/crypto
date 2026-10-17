#!/bin/sh
# Thin orchestration for PROTOCOL.md (no computation here).
#
#   run_arms.sh <raw_dir> <baseline_bin_dir> <candidate_bin_dir> <part>
#
# <part> is one of: bsgs1 (BSGS, 1 thread), bsgs3 (BSGS, 3 threads, 48 bits),
# rho1 (rho, 1 thread), rho3 (rho, 3 threads).  Every timed process runs under
# tools/isolated_bench.py; its isolation record goes to <raw_dir>/isolation.jsonl
# and the binary's own JSON line, tagged by --label, to <raw_dir>/<part>.jsonl.
# Arms per round are rotated A B A2 / B A2 A / A2 A B; A and A2 are the
# baseline binary (the A/A control), B the candidate.
set -eu
raw=$1; base=$2; cand=$3; part=$4
here=$(cd "$(dirname "$0")" && pwd)
study=$(cd "$here/.." && pwd)
repo=$(cd "$study/../.." && pwd)
curves="$study/inputs/curves.json"
iso="$raw/isolation.jsonl"
out="$raw/$part.jsonl"
mkdir -p "$raw"

# run <cpus> <label> cmd...: retry while the preflight finds the host busy,
# wait for the lock while a build holds it (isolated_bench.py busy), and skip
# a label already recorded, so an interrupted part resumes where it stopped.
run() {
  cpus=$1; label=$2; shift 2
  if [ -f "$out" ] && grep -q "\"label\":\"$label\"" "$out"; then
    return 0
  fi
  tries=0
  while :; do
    if python3 "$repo/tools/isolated_bench.py" run --wait --cpus "$cpus" --out "$iso" \
        --label "$label" -- "$@" > "$raw/.last" 2> "$raw/.err"; then
      cat "$raw/.last" >> "$out"
      return 0
    fi
    if grep -q "machine is busy" "$raw/.err" && [ "$tries" -lt 60 ]; then
      tries=$((tries + 1)); sleep 15; continue
    fi
    cat "$raw/.err" >&2
    return 1
  done
}

arms_for_round() {
  case $(( ($1 - 1) % 3 )) in
    0) echo "A B A2" ;; 1) echo "B A2 A" ;; *) echo "A2 A B" ;;
  esac
}

bin_for() { # bin_for <arm> <name>
  case $1 in B) echo "$cand/$2" ;; *) echo "$base/$2" ;; esac
}

rounds_for() { if [ "$1" -ge 52 ]; then echo 5; else echo 7; fi; }

case $part in
  bsgs1|bsgs3)
    if [ "$part" = bsgs1 ]; then threads=1; cpus=3; sizes="40 44 48 52"; else threads=3; cpus=1-3; sizes=48; fi
    for bits in $sizes; do
      c=$(jq -c ".curves[] | select(.role == \"bsgs\" and .bits == $bits)" "$curves")
      p=$(echo "$c" | jq -r .p); a=$(echo "$c" | jq -r .a); b=$(echo "$c" | jq -r .b)
      gx=$(echo "$c" | jq -r .generator.x); gy=$(echo "$c" | jq -r .generator.y)
      n=$(echo "$c" | jq -r .subgroup_order)
      r=1; R=$(rounds_for "$bits")
      while [ "$r" -le "$R" ]; do
        for kern in chain paired; do
          for arm in $(arms_for_round "$r"); do
            label="$part-$bits-$kern-$arm-r$r"
            run "$cpus" "$label" env RAYON_NUM_THREADS=$threads BSGS_STEP=$kern \
              "$(bin_for "$arm" lit2026_w11_bsgs_time)" "$p" "$a" "$b" "$gx" "$gy" "$n" 1 "$bits" "$label"
          done
        done
        r=$((r + 1))
      done
    done
    ;;
  rho1|rho3)
    if [ "$part" = rho1 ]; then threads=1; cpus=3; else threads=3; cpus=1-3; fi
    for bits in 32 36 40 44 48 52; do
      r=1; R=$(rounds_for "$bits")
      while [ "$r" -le "$R" ]; do
        for arm in $(arms_for_round "$r"); do
          label="$part-$bits-$arm-r$r"
          run "$cpus" "$label" env RAYON_NUM_THREADS=$threads \
            "$(bin_for "$arm" generic_dlp_rho_time)" "$bits" "$threads" 1 "$label"
        done
        r=$((r + 1))
      done
    done
    ;;
  *) echo "unknown part $part" >&2; exit 2 ;;
esac
rm -f "$raw/.last" "$raw/.err"

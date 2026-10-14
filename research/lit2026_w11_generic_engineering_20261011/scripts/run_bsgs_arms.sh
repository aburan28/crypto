#!/bin/sh
# lit2026 W11 (a): interleaved BSGS timing arms (thin orchestration only).
#
#   run_bsgs_arms.sh <out.jsonl> <rounds> <orig_bin> <cand_bin> <pair_k> <curves.jsonl>
#
# curves.jsonl: one line per curve as printed by
#   lit2026_w11_challenge_ladder derive  (fields bits, p, n, Qx, Qy).
# Arms per round, rotated: original (baseline binary), previous
# (BSGS_STEP=chain), paired, coz.  Every sample line is the binary's JSON
# with the arm, the round and `uptime` appended.
set -eu
out=$1; rounds=$2; orig=$3; cand=$4; k=$5; curves=$6
: > "$out"
targets_for() {
  case $1 in
    40) echo 64 ;; 44) echo 32 ;; 48) echo 16 ;; 52) echo 8 ;; *) echo 4 ;;
  esac
}
arms="original previous paired coz"
r=1
while [ "$r" -le "$rounds" ]; do
  while read -r line; do
    bits=$(printf '%s' "$line" | sed -E 's/.*"bits":([0-9]+).*/\1/')
    p=$(printf '%s' "$line" | sed -E 's/.*"p":([0-9]+).*/\1/')
    n=$(printf '%s' "$line" | sed -E 's/.*"n":([0-9]+).*/\1/')
    gx=$(printf '%s' "$line" | sed -E 's/.*"Qx":([0-9]+).*/\1/')
    gy=$(printf '%s' "$line" | sed -E 's/.*"Qy":([0-9]+).*/\1/')
    t=$(targets_for "$bits")
    # rotate the arm order by round
    set -- $arms
    i=0
    while [ $i -lt $(( (r - 1) % 4 )) ]; do
      first=$1; shift; set -- "$@" "$first"; i=$((i + 1))
    done
    for arm in "$@"; do
      up=$(uptime | sed -E 's/.*load averages?: //')
      case $arm in
        original) res=$("$orig" "$p" 0 7 "$gx" "$gy" "$n" "$t" "$bits" "secp${bits}k1") ;;
        previous) res=$(BSGS_STEP=chain "$cand" "$p" 0 7 "$gx" "$gy" "$n" "$t" "$bits" "secp${bits}k1") ;;
        paired) res=$(BSGS_STEP=paired BSGS_PAIR_K=$k "$cand" "$p" 0 7 "$gx" "$gy" "$n" "$t" "$bits" "secp${bits}k1") ;;
        coz) res=$(BSGS_STEP=coz BSGS_PAIR_K=$k "$cand" "$p" 0 7 "$gx" "$gy" "$n" "$t" "$bits" "secp${bits}k1") ;;
      esac
      printf '%s\n' "$res" | sed -E "s/}\$/,\"arm\":\"$arm\",\"round\":$r,\"uptime_load\":\"$up\"}/" >> "$out"
    done
  done < "$curves"
  r=$((r + 1))
done

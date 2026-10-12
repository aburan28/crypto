#!/bin/sh
# Five identical-command A/A pairs for both native arms; run from repo root.
set -eu
note=research/notes/ecc2k130/n37_native_one_target_20261003
points=research/notes/ecc2k130/disjoint_cold_v2_20261001/fixtures/n37_L1024_b02.points.jsonl
public_point=$(sed -n '1p' "$points" | jq -r 'join(",")')
mkdir -p "$note/aa"
j=0
while [ "$j" -lt 5 ]; do
    for copy in a b; do
        prefix="$note/aa/q00_r${j}_${copy}"
        if [ -e "$prefix.replay.json" ]; then
            continue
        fi
        if [ -e "$prefix.ic.json" ] || [ -e "$prefix.rho.jsonl" ]; then
            echo "partial A/A artifact at $prefix; preserve and inspect it" >&2
            exit 1
        fi
        /usr/bin/time -p env RAYON_NUM_THREADS=1 \
            target/release/examples/n37_native_m6_one_target 0 "$prefix.ic.json" \
            >"$prefix.ic.stdout" 2>"$prefix.ic.stderr"
        /usr/bin/time -p env RAYON_NUM_THREADS=1 KIC_RHO_RUNG=3 \
            KIC_RHO_LANES=32 KIC_RHO_DP_BITS=8 \
            KIC_RHO_BATCH_CORPUS=compact-disjoint-cold-v2-n37-L1024-b02-20261001 \
            KIC_RHO_TARGET_POINT="$public_point" \
            target/release/examples/koblitz_rho_batch_ks_strong_online \
            37 0 signed_frobenius 1 2026100110102 \
            >"$prefix.rho.jsonl" 2>"$prefix.rho.stderr"
        target/release/examples/n37_native_m6_one_target_replay \
            "$prefix.ic.json" "$prefix.rho.jsonl" "$prefix.replay.json" \
            >"$prefix.replay.stdout" 2>"$prefix.replay.stderr"
        echo "q00 A/A round $j copy $copy complete"
    done
    j=$((j + 1))
done

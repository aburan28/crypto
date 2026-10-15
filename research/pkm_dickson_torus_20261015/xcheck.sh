#!/bin/bash
# Build-provenance cross-check (README §1): re-measure a subset of the
# recorded cells with the binary rebuilt in this worktree's own target dir,
# same flags, and compare every deterministic field with the recorded rows
# (pkm_tower_check compare-rows; only ms and wall_s are left out).
set -u
here=$(cd "$(dirname "$0")" && pwd)
bin=${BIN:?pkm_tower_pilot}; chk=${CHK:?pkm_tower_check}
cd "$here/xcheck"
ulimit -v 2621440
export RAYON_NUM_THREADS=1
common=(--engine tower --budget 7200 --max-nnz 2000000000 --ladder-t none --planted 0 --random 2)
x() { name=$1; ref=$2; shift 2
  nice -n 10 "$bin" "${common[@]}" "$@" --out "$name.jsonl" 2> "$name.log"
  echo "== $name vs $ref (exit $?)"; "$chk" compare-rows "../runs/$ref.jsonl" "$name.jsonl" 2>&1 | tail -2; }
x X0 R0-isogeny-m2-p0-poscontrol --p 786433 --kinds isogeny --m 2 --controls tower --t-min 1 --max-t 5
x X1 R1-torus-m2-T31 --p 2046820351 --kinds dickson-torus --m 2 --controls tower --t-min 1 --max-t 9
x X2 R2-torus-m2-T24 --p 14680063 --kinds dickson-torus --m 2 --controls tower --t-min 1 --max-t 8
x X3 R3-split-dickson-m2-p1 --p 2013265921 --kinds dickson --m 2 --controls tower --t-min 1 --max-t 8
x X4 R4-torus-m2-T31-null --p 2046820351 --kinds dickson-torus --m 2 --controls null --t-min 1 --max-t 8
x X6 R6-torus-m34-T31 --p 2046820351 --kinds dickson-torus --m 3,4 --t-min 2 --max-t-m3 4 --max-t-m4 3 --controls tower
x X8 R8-split-dickson-m34-p1 --p 2013265921 --kinds dickson --m 3,4 --t-min 2 --max-t-m3 4 --max-t-m4 3 --controls tower
x X10 R10-torus-m34-T24-null --p 14680063 --kinds dickson-torus --m 3,4 --t-min 2 --max-t-m3 4 --max-t-m4 3 --controls null --null-chain
x X11 R11-kummer-p2-deep --p 3221225473 --kinds kummer --m 2,3 --controls tower --t-min 1 --max-t 8 --max-t-m3 4
x X11n R11n-kummer-m2-p2-deep-null --p 3221225473 --kinds kummer --m 2 --controls null --t-min 1 --max-t 8

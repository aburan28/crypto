#!/usr/bin/env bash
# Separate precomputed-tree retry after all table campaigns and the batched run.
set -u
study=research/conductor_gap_cm_orbits_20261015
binary=/workspace/native-deps/cm-sqrt-velu-precomputed-final
sha256sum "$binary" "$study/native/sqrt_velu.cpp" \
  "$study/native/run_96_precomputed.sh" "$study/prepared/large-batched-1007.txt" \
  /workspace/native-deps/root/usr/lib/x86_64-linux-gnu/libflint.so.19 \
  /workspace/native-deps/root/usr/bin/time > "$study/results/precomputed-96-pre-run.sha256"
if [ $? -ne 0 ]; then exit 2; fi
/workspace/native-deps/root/usr/bin/time -v -o "$study/results/precomputed-96-time.txt" \
  timeout --signal=TERM --kill-after=10s 3600s \
  prlimit --as=15032385536 -- "$binary" \
  < "$study/prepared/large-batched-1007.txt" \
  > "$study/results/precomputed-96.jsonl" 2> "$study/results/precomputed-96.stderr"
status=$?
printf '{"case":1007,"exit":%s,"wall_cap_seconds":3600,"address_space_cap_bytes":15032385536,"isolation":"L0; no concurrent solver tables","scope":"same exact input; precomputed remainder tree, single giant-point batch"}\n' "$status" > "$study/results/precomputed-96-exit.json"

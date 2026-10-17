#!/usr/bin/env bash
# Separate optimized, bounded-batch retry. Refuse concurrent local compilers.
set -u
study=research/conductor_gap_cm_orbits_20261015
binary=/workspace/native-deps/cm-sqrt-velu-v4
ps -eo pid=,stat=,comm= > "$study/results/v4-96-prelaunch-processes.txt"
if awk '$2 !~ /^Z/ && $3 ~ /^(rustc|clippy-driver|cargo|cargo-clippy|cc1plus|ecbench|conductor_gap_.*|cm-sqrt.*|cm-folded.*)$/ { found=1 } END { exit !found }' "$study/results/v4-96-prelaunch-processes.txt"; then
  echo 'Refusing construction while local compiler or solver remains active' >&2
  exit 75
fi
cat /sys/fs/cgroup/memory.events > "$study/results/v4-96-memory-events-before.txt"
sha256sum "$binary" "$study/native/sqrt_velu_batched_v4.cpp" \
  "$study/native/sqrt_velu.cpp" "$study/native/run_96_v4_solo.sh" \
  "$study/prepared/large-batched-1007.txt" \
  /workspace/native-deps/root/usr/lib/x86_64-linux-gnu/libflint.so.19 \
  /workspace/native-deps/root/usr/bin/time > "$study/results/v4-96-pre-run.sha256"
if [ $? -ne 0 ]; then exit 2; fi
/workspace/native-deps/root/usr/bin/time -v -o "$study/results/v4-96-time.txt" \
  timeout --signal=TERM --kill-after=10s 5400s \
  prlimit --as=11811160064 -- "$binary" 2097152 \
  < "$study/prepared/large-batched-1007.txt" \
  > "$study/results/v4-96.jsonl" 2> "$study/results/v4-96.stderr"
status=$?
cat /sys/fs/cgroup/memory.events > "$study/results/v4-96-memory-events-after.txt"
printf '{"case":1007,"exit":%s,"wall_cap_seconds":5400,"address_space_cap_bytes":11811160064,"batch_cap":2097152,"isolation":"L0; no concurrent local compilers or solver tables","scope":"same exact input; closed-form moment coefficients and delayed reductions"}\n' "$status" > "$study/results/v4-96-exit.json"

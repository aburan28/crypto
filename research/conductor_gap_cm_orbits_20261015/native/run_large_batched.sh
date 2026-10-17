#!/usr/bin/env bash
# Bounded orchestration; every case, failure and censor gets a receipt.
set -u
study=research/conductor_gap_cm_orbits_20261015
test -x /workspace/native-deps/root/usr/bin/time || exit 127
for case_id in 1005 1006 1007; do
  awk -v key="$case_id" '$1==key' "$study/prepared/new-kernels.txt" > "$study/prepared/large-batched-$case_id.txt"
done
sha256sum /workspace/native-deps/cm-sqrt-velu-v3b /workspace/native-deps/root/usr/bin/time "$study/native/run_large_batched.sh" "$study/native/sqrt_velu_batched.cpp" "$study/native/sqrt_velu.cpp" "$study"/prepared/large-batched-*.txt > "$study/results/batched-large-pre-run.sha256"
for case_id in 1005 1006 1007; do
  admitted=0
  for attempt in $(seq 1 40); do
    anon_bytes=$(awk '$1=="anon" {print $2}' /sys/fs/cgroup/memory.stat)
    if [ "$anon_bytes" -lt 10000000000 ]; then admitted=1; break; fi
    sleep 30
  done
  if [ "$admitted" -eq 0 ]; then
    printf '{"case":%s,"status":"not-admitted-memory","wall_cap_seconds":1800,"address_space_cap_bytes":6442450944}\n' "$case_id" > "$study/results/large-batched-$case_id-exit.json"
    continue
  fi
  /workspace/native-deps/root/usr/bin/time -v -o "$study/results/large-batched-$case_id-time.txt" \
    timeout --signal=TERM --kill-after=10s 1800s \
    prlimit --as=6442450944 -- /workspace/native-deps/cm-sqrt-velu-v3b 1048576 \
    < "$study/prepared/large-batched-$case_id.txt" \
    > "$study/results/large-batched-$case_id.jsonl" \
    2> "$study/results/large-batched-$case_id.stderr"
  status=$?
  printf '{"case":%s,"exit":%s,"wall_cap_seconds":1800,"address_space_cap_bytes":6442450944,"batch_cap":1048576,"isolation":"L0"}\n' "$case_id" "$status" > "$study/results/large-batched-$case_id-exit.json"
done

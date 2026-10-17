#!/usr/bin/env bash
# Separate retry after the co-scheduled 6 GiB campaign. Same exact 96-bit input.
set -u
study=research/conductor_gap_cm_orbits_20261015
sha256sum /workspace/native-deps/cm-sqrt-velu-v3b /workspace/native-deps/root/usr/bin/time \
  "$study/native/run_96_expanded.sh" "$study/native/sqrt_velu_batched.cpp" \
  "$study/native/sqrt_velu.cpp" "$study/prepared/large-batched-1007.txt" \
  /workspace/native-deps/root/usr/lib/x86_64-linux-gnu/libflint.so.19 \
  > "$study/results/expanded-96-pre-run.sha256"
if [ $? -ne 0 ]; then exit 2; fi
/workspace/native-deps/root/usr/bin/time -v -o "$study/results/expanded-96-time.txt" \
  timeout --signal=TERM --kill-after=10s 3600s \
  prlimit --as=11811160064 -- /workspace/native-deps/cm-sqrt-velu-v3b 1048576 \
  < "$study/prepared/large-batched-1007.txt" \
  > "$study/results/expanded-96.jsonl" 2> "$study/results/expanded-96.stderr"
status=$?
printf '{"case":1007,"exit":%s,"wall_cap_seconds":3600,"address_space_cap_bytes":11811160064,"batch_cap":1048576,"isolation":"L0; validation may overlap","scope":"same exact input as earlier 6 GiB attempt"}\n' "$status" > "$study/results/expanded-96-exit.json"

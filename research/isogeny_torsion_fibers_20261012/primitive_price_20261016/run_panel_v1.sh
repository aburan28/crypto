#!/usr/bin/env bash
# Run from the repository root. Frozen protocol: PROTOCOL_V1.md.
set -u

out="${1:?usage: run_panel_v1.sh NEW_OUTPUT_DIRECTORY}"
test ! -e "$out" || { echo "output already exists: $out" >&2; exit 97; }
base=research/isogeny_torsion_fibers_20261012/portable_order_trend_20261015/order1000
public="$base/BATCH_1237_1201_V2.public.json"
answers="$base/BATCH_1237_1201_V2.checker.json"
binary=target/release/examples/isogeny_fiber_primitive_price
checker=target/release/examples/different_j_inverse_velu_check
meter=target/release/examples/different_j_order_process_meter

test "$(shasum -a 256 "$public" | awk '{print $1}')" = \
  24f2c56a9ba56f91c954f38b7752edf24fee5ed027b49b1ac0747771945099c4 || exit 99
mkdir -p "$out/raw"
git rev-parse HEAD > "$out/git-head.txt"
date -u '+%Y-%m-%dT%H:%M:%SZ' > "$out/started-utc.txt"
uname -a > "$out/host.txt"
rustc --version >> "$out/host.txt"
shasum -a 256 "$public" "$answers" "$binary" "$checker" "$meter" \
  examples/isogeny_fiber_primitive_price.rs \
  research/isogeny_torsion_fibers_20261012/primitive_price_20261016/run_panel_v1.sh \
  > "$out/build-sha256.txt"
: > "$out/run.jsonl"
: > "$out/meter.jsonl"
arms=(inverse-velu full-image source-bsgs rho-negation rho-negation-aa)
any_failed=0
for prefix in 1 4 16 64; do
  for round in 0 1 2; do
    for offset in 0 1 2 3 4; do
      arm=${arms[$(((round + offset) % 5))]}
      child="$out/raw/prefix-$prefix-round-$round-$arm"
      gtimeout 300 "$meter" "$binary" "$public" "$round" "$arm" "$prefix" \
        "$child.stdout" "$child.stderr" "$child.meter.json" \
        > "$child.meter.stdout" 2> "$child.meter.stderr"
      child_exit=$?
      printf '%s\n' "$child_exit" > "$child.exit"
      if (( child_exit != 0 )); then any_failed=1; fi
      if [[ -e "$child.stdout" ]]; then cat "$child.stdout" >> "$out/run.jsonl"; fi
      if [[ -e "$child.meter.json" ]]; then
        cat "$child.meter.json" >> "$out/meter.jsonl"
        printf '\n' >> "$out/meter.jsonl"
      fi
    done
  done
done
if (( any_failed == 0 )); then
  gtimeout 300 "$checker" "$public" "$answers" "$out/run.jsonl" "$out/meter.jsonl" \
    > "$out/replay.json" 2> "$out/replay.stderr"
  replay_exit=$?
  printf '%s\n' "$replay_exit" > "$out/replay.exit"
  if (( replay_exit != 0 )); then any_failed=1; fi
else
  printf 'skipped_after_failed_children\n' > "$out/replay.status"
fi
printf '%s\n' "$any_failed" > "$out/exit"
(cd "$out" && find . -type f ! -name ARTIFACT_FILES.sha256 -print0 | sort -z | xargs -0 shasum -a 256) \
  > "$out/ARTIFACT_FILES.sha256"
exit "$any_failed"

#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
cd "$repo"
out=research/f6_ic_rho_online_gap_20261006
source=research/f6_ic_compact_refuted_20261006
archived=research/f6_ic_geometric_closure_20261003/eight_target
binary=target/bench_bins/f6_ic_compact_refuted_worker
mkdir -p "$out/inputs" "$out/candidates" "$out/workloads"
sha() { shasum -a 256 "$1" | awk '{print $1}'; }
test "$(sha "$binary")" = "$(sed -n 's/^binary_sha256=//p' "$source/BUILD_IDENTITY.txt")"
shasum -a 256 -c "$source/SOURCE_SHA256SUMS"
cp "$source/candidates/compact_baseline.json" "$out/candidates/f6_ic.json"
cp "$source/BUILD_IDENTITY.txt" "$out/BUILD_IDENTITY.txt"
cp "$source/SOURCE_SHA256SUMS" "$out/SOURCE_SHA256SUMS"
jq -S -c --arg bin "$(sha "$binary")" \
  --arg worker "$(sha examples/ic_tournament_worker.rs)" \
  --arg engine "$(sha src/cryptanalysis/koblitz_index_calculus.rs)" \
  '{schema_version:1,algorithm:"signed-frobenius-batched-affine",
    curve_id:"EC1N17Ckb1hbbe2b5b6b1e6",binary_sha256:$bin,
    worker_source_sha256:$worker,engine_source_sha256:$engine,
    seed:20261004039,requested_walks:1,jump_count:16,
    max_restarts:64,max_iterations_per_restart:65536,
    progress_interval:256,walk_and_collision_policy:"source-default-v1",
    target_input:"supplied_public_point",online_boundary:
    "after-reusable-fast-curve-preparation-before-target-lift"}' \
  -n > "$out/rho_reference.json"
printf 'target\tworkload_id\tic_candidate_id\tic_input_sha256\trho_input_sha256\tpoint_sha256\n' > "$out/FREEZE.tsv"
for target in 1 7; do
    ic="$out/inputs/T$target-ic.json"
    rho="$out/inputs/T$target-rho.json"
    cp "$source/inputs/T$target-compact_baseline.json" "$ic"
    cp "$archived/T$target/workload.json" "$out/workloads/T$target.json"
    jq -S -c '.mode="rho" | del(.prepared) |
        .config.solver="rho" | .config.max_trials=65536 |
        .config.rho_parallel_walks=1' "$ic" > "$rho"
    point_sha=$(jq -S -c '.public_targets[0]' "$ic" | shasum -a 256 | cut -d ' ' -f 1)
    test "$(jq -S -c '.public_targets[0]' "$ic")" = "$(jq -S -c '.public_targets[0]' "$rho")"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$target" \
      "$(jq -r .workload_id "$out/workloads/T$target.json")" \
      "$(jq -r .candidate_id "$out/candidates/f6_ic.json")" \
      "$(sha "$ic")" "$(sha "$rho")" "$point_sha" >> "$out/FREEZE.tsv"
done
shasum -a 256 "$out/rho_reference.json" "$out/candidates/f6_ic.json" \
  "$out/workloads/T1.json" "$out/workloads/T7.json" > "$out/FROZEN_SHA256SUMS"

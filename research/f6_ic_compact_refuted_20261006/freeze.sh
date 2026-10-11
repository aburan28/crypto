#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
cd "$repo"
out=research/f6_ic_compact_refuted_20261006
base=research/f6_ic_inherited_basis_profile_20261006
archived=research/f6_ic_geometric_closure_20261003/eight_target
mkdir -p "$out/inputs" "$out/candidates" target/bench_bins
cp target/release/examples/ic_tournament_worker target/bench_bins/f6_ic_compact_refuted_worker
binary=target/bench_bins/f6_ic_compact_refuted_worker
sha() { shasum -a 256 "$1" | awk '{print $1}'; }
commit=$(git rev-parse HEAD)
worker_hash=$(sha examples/ic_tournament_worker.rs)
groebner_hash=$(sha src/cryptanalysis/koblitz_groebner.rs)
inherited_hash=$(sha src/cryptanalysis/inherited_f4.rs)
engine_hash=$(sha src/cryptanalysis/koblitz_index_calculus.rs)
measurement_hash=$(sha src/cryptanalysis/ic_measurement.rs)
binary_hash=$(sha "$binary")
source_manifest="$base/candidates/f6_ic_basis_profile.json"
old_worker=$(jq -r '.record.implementation.components[]|select(.role=="ic-worker")|.sha256' "$source_manifest")
old_groebner=$(jq -r '.record.implementation.components[]|select(.role=="boolean-groebner")|.sha256' "$source_manifest")
old_engine=$(jq -r '.record.implementation.components[]|select(.role=="ic-engine")|.sha256' "$source_manifest")
old_measurement=$(jq -r '.record.implementation.components[]|select(.role=="exclusive-measurement")|.sha256' "$source_manifest")
printf 'commit=%s\nbinary_sha256=%s\n' "$commit" "$binary_hash" > "$out/BUILD_IDENTITY.txt"
shasum -a 256 examples/ic_tournament_worker.rs src/cryptanalysis/koblitz_groebner.rs \
  src/cryptanalysis/inherited_f4.rs \
  src/cryptanalysis/koblitz_index_calculus.rs src/cryptanalysis/ic_measurement.rs \
  research/ic_candidate_tournament_20260915/goal_20260924/prepared-f5-runtime-v1/mathematics.json \
  > "$out/SOURCE_SHA256SUMS"
for arm in compact_baseline compact_candidate; do
    source="$base/candidates/f6_ic_basis_profile.json"
    draft="$out/candidates/$arm.draft.json"
    jq -S -c \
      --arg commit "$commit" --arg binary "$binary_hash" \
      --arg ow "$old_worker" --arg nw "$worker_hash" \
      --arg og "$old_groebner" --arg ng "$groebner_hash" \
      --arg inherited "$inherited_hash" \
      --arg oe "$old_engine" --arg ne "$engine_hash" \
      --arg om "$old_measurement" --arg nm "$measurement_hash" \
      --arg arm "$arm" \
      'walk(if type=="string" then
          if .==$ow then $nw elif .==$og then $ng elif .==$oe then $ne
          elif .==$om then $nm else . end else . end)
       | .record.implementation.code_snapshot_commit=$commit
       | .record.implementation.components += [{"role":"inherited-f4","sha256":$inherited}]
       | .record.implementation.flags.binary_sha256=$binary
       | .record.implementation.flags.profiling="none"
       | .record.implementation.flags.active_multipliers=false
       | .record.implementation.flags.direct_fused_pack=false
       | .record.implementation.flags.support_local_stream=false
       | .record.implementation.flags.support_local_bitmap_columns=false
       | .record.implementation.flags.support_local_profile=false
       | .record.implementation.flags.inherited_basis_profile=false
       | .record.implementation.flags.compact_refuted=($arm=="compact_candidate")
       | .record.point_decomposition.cache_policy="process-local-template-and-layout-caches-source-defaults"
       | .candidate_id=null' "$source" > "$draft"
    digest=$(jq -S -c '.record' "$draft" | tr -d '\n' | shasum -a 256 | cut -c1-12)
    prefix=$(jq -r .candidate_id "$source" | sed -E 's/h[0-9a-f]{12}$//')
    candidate_id="${prefix}h${digest}"
    jq -S -c --arg id "$candidate_id" '.candidate_id=$id' "$draft" > "$out/candidates/$arm.json"
    rm "$draft"
done
printf 'target\tarm\tinput_sha256\tcandidate_id\tcandidate_sha256\n' > "$out/FREEZE.tsv"
for target in 1 7; do
    for arm in compact_baseline compact_candidate; do
        jq -S -c \
          --arg arm "$arm" \
          '.config.active_multipliers=false | .config.direct_fused_pack=false
           | .config.support_local_stream=false
           | .config.support_local_bitmap_columns=false
           | .config.support_local_profile=false
           | .config.inherited_basis_profile=false
           | .config.compact_refuted=($arm=="compact_candidate")' \
          "$archived/T$target/input-f6_ic.json" > "$out/inputs/T$target-$arm.json"
        input="$out/inputs/T$target-$arm.json"
        candidate="$out/candidates/$arm.json"
        printf '%s\t%s\t%s\t%s\t%s\n' "$target" "$arm" "$(sha "$input")" \
          "$(jq -r .candidate_id "$candidate")" "$(sha "$candidate")" >> "$out/FREEZE.tsv"
    done
done

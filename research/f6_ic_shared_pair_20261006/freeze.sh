#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
cd "$repo"
out=research/f6_ic_shared_pair_20261006
base=research/f6_ic_adaptive_pair_20261005
archived=research/f6_ic_geometric_closure_20261003/eight_target
mkdir -p "$out/inputs" "$out/candidates" target/bench_bins
cp target/release/examples/ic_tournament_worker target/bench_bins/f6_ic_shared_pair_worker
binary=target/bench_bins/f6_ic_shared_pair_worker
sha() { shasum -a 256 "$1" | awk '{print $1}'; }
commit=$(git rev-parse HEAD)
worker_hash=$(sha examples/ic_tournament_worker.rs)
groebner_hash=$(sha src/cryptanalysis/koblitz_groebner.rs)
engine_hash=$(sha src/cryptanalysis/koblitz_index_calculus.rs)
measurement_hash=$(sha src/cryptanalysis/ic_measurement.rs)
binary_hash=$(sha "$binary")
source_manifest="$base/candidates/f6_ic_pair.json"
old_worker=$(jq -r '.record.implementation.components[]|select(.role=="ic-worker")|.sha256' "$source_manifest")
old_groebner=$(jq -r '.record.implementation.components[]|select(.role=="boolean-groebner")|.sha256' "$source_manifest")
old_engine=$(jq -r '.record.implementation.components[]|select(.role=="ic-engine")|.sha256' "$source_manifest")
old_measurement=$(jq -r '.record.implementation.components[]|select(.role=="exclusive-measurement")|.sha256' "$source_manifest")
printf 'commit=%s\nbinary_sha256=%s\n' "$commit" "$binary_hash" > "$out/BUILD_IDENTITY.txt"
shasum -a 256 examples/ic_tournament_worker.rs src/cryptanalysis/koblitz_groebner.rs \
  src/cryptanalysis/koblitz_index_calculus.rs src/cryptanalysis/ic_measurement.rs \
  research/ic_candidate_tournament_20260915/goal_20260924/prepared-f5-runtime-v1/mathematics.json \
  > "$out/SOURCE_SHA256SUMS"
for arm in inherited_f4 f6_ic f6_ic_pair f6_ic_shared_pair; do
    source_arm=$arm
    [ "$arm" = f6_ic_shared_pair ] && source_arm=f6_ic_pair
    source="$base/candidates/$source_arm.json"
    draft="$out/candidates/$arm.draft.json"
    jq -S -c \
      --arg commit "$commit" --arg binary "$binary_hash" \
      --arg ow "$old_worker" --arg nw "$worker_hash" \
      --arg og "$old_groebner" --arg ng "$groebner_hash" \
      --arg oe "$old_engine" --arg ne "$engine_hash" \
      --arg om "$old_measurement" --arg nm "$measurement_hash" \
      --arg arm "$arm" \
      'walk(if type=="string" then
          if .==$ow then $nw elif .==$og then $ng elif .==$oe then $ne
          elif .==$om then $nm else . end else . end)
       | .record.implementation.code_snapshot_commit=$commit
       | .record.implementation.flags.binary_sha256=$binary
       | .record.implementation.flags.f6_pair_index=($arm=="f6_ic_pair")
       | .record.implementation.flags.f6_shared_pair_index=($arm=="f6_ic_shared_pair")
       | if $arm=="f6_ic_shared_pair" then
           .record.point_decomposition.solver="f6_ic_shared_pair"
           | .record.point_decomposition.geometric_oracle.version=5
           | .record.point_decomposition.geometric_oracle.closure="m3-one-fixed-x-adaptive-exact-unordered-pair-sum-shared-on-base"
           | .record.point_decomposition.geometric_oracle.cache_policy="per-process-factor-base-and-curve-fingerprint;first-build-charged-online"
         else . end
       | .candidate_id=null' "$source" > "$draft"
    digest=$(jq -S -c '.record' "$draft" | tr -d '\n' | shasum -a 256 | cut -c1-12)
    prefix=$(jq -r .candidate_id "$source" | sed -E 's/h[0-9a-f]{12}$//')
    candidate_id="${prefix}h${digest}"
    jq -S -c --arg id "$candidate_id" '.candidate_id=$id' "$draft" > "$out/candidates/$arm.json"
    rm "$draft"
done
printf 'target\tarm\tinput_sha256\tcandidate_id\tcandidate_sha256\n' > "$out/FREEZE.tsv"
for target in 1 7; do
    for arm in inherited_f4 f6_ic f6_ic_pair f6_ic_shared_pair; do
        if [ "$arm" = f6_ic_shared_pair ]; then
            jq -S -c '.config.solver="f6_ic_shared_pair"' \
              "$base/inputs/T$target-f6_ic_pair.json" > "$out/inputs/T${target}-$arm.json"
            input="$out/inputs/T${target}-$arm.json"
        elif [ "$arm" = f6_ic_pair ]; then
            input="$base/inputs/T$target-$arm.json"
        else
            input="$archived/T$target/input-$arm.json"
        fi
        candidate="$out/candidates/$arm.json"
        printf '%s\t%s\t%s\t%s\t%s\n' "$target" "$arm" "$(sha "$input")" \
          "$(jq -r .candidate_id "$candidate")" "$(sha "$candidate")" >> "$out/FREEZE.tsv"
    done
done

#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
cd "$repo"
out=research/f6_ic_specialise_breakdown_20261006
prior=research/f6_ic_inherited_basis_profile_20261006
binary=target/bench_bins/f6_ic_specialise_breakdown_worker
mkdir -p "$out/inputs" "$out/candidates" target/bench_bins
cp target/release/examples/ic_tournament_worker "$binary"
sha() { shasum -a 256 "$1" | awk '{print $1}'; }
printf 'commit=%s\nbinary_sha256=%s\n' "$(git rev-parse HEAD)" "$(sha "$binary")" > "$out/BUILD_IDENTITY.txt"
shasum -a 256 examples/ic_tournament_worker.rs src/cryptanalysis/inherited_f4.rs \
    src/cryptanalysis/koblitz_groebner.rs \
    src/cryptanalysis/koblitz_index_calculus.rs \
    src/cryptanalysis/ic_measurement.rs \
    research/ic_candidate_tournament_20260915/goal_20260924/prepared-f5-runtime-v1/mathematics.json \
    > "$out/SOURCE_SHA256SUMS"
worker_sha=$(sha examples/ic_tournament_worker.rs)
inherited_sha=$(sha src/cryptanalysis/inherited_f4.rs)
binary_sha=$(sha "$binary")
commit=$(git rev-parse HEAD)
prefix=IC1N17Ckb1fb62PDP3f6RCsampleLAgaussTDpdpISO0h
printf 'target\tarm\tinput_sha256\tcandidate_id\tcandidate_sha256\n' > "$out/FREEZE.tsv"
for target in 1 7; do
    input="$out/inputs/T${target}-f6_ic_specialise_breakdown.json"
    candidate="$out/candidates/T${target}-f6_ic_specialise_breakdown.json"
    jq '.config.specialise_profile = true' "$prior/inputs/T${target}-f6_ic_basis_profile.json" > "$input"
    record="$out/candidates/T${target}.record.json"
    jq --arg commit "$commit" --arg worker_sha "$worker_sha" \
        --arg inherited_sha "$inherited_sha" --arg binary_sha "$binary_sha" \
        '.record | .implementation.code_snapshot_commit = $commit |
         .implementation.components |= ([.[] | if .role == "ic-worker" then .sha256 = $worker_sha else . end] + [{"role":"inherited-f4","sha256":$inherited_sha}]) |
         .implementation.flags.binary_sha256 = $binary_sha |
         .implementation.flags.specialise_profile = true |
         .implementation.flags.profiling = "target-online-inherited-specialise-breakdown-ns-v1"' \
         "$prior/candidates/f6_ic_basis_profile.json" > "$record"
    digest=$(jq -S -c . "$record" | tr -d '\n' | shasum -a 256 | cut -c1-12)
    candidate_id="${prefix}${digest}"
    jq -n --arg id "$candidate_id" --slurpfile record "$record" \
        '{candidate_id:$id,record:$record[0]}' > "$candidate"
    rm "$record"
    printf '%s\tf6_ic_specialise_breakdown\t%s\t%s\t%s\n' \
        "$target" "$(sha "$input")" "$candidate_id" "$(sha "$candidate")" >> "$out/FREEZE.tsv"
done

#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
cd "$repo"
out=research/f6_ic_specialise_breakdown_20261006
archived=research/f6_ic_geometric_closure_20261003/eight_target
binary=target/bench_bins/f6_ic_specialise_breakdown_worker
sha() { shasum -a 256 "$1" | awk '{print $1}'; }
expected_binary=$(sed -n 's/^binary_sha256=//p' "$out/BUILD_IDENTITY.txt")
test "$(sha "$binary")" = "$expected_binary"
shasum -a 256 -c "$out/SOURCE_SHA256SUMS"
while IFS="$(printf '\t')" read -r target arm input_sha candidate_id candidate_sha; do
    [ "$target" = target ] && continue
    input="$out/inputs/T${target}-$arm.json"
    candidate="$out/candidates/T${target}-$arm.json"
    test "$(sha "$input")" = "$input_sha"
    test "$(sha "$candidate")" = "$candidate_sha"
    test "$(jq -r .candidate_id "$candidate")" = "$candidate_id"
    digest=$(jq -S -c '.record' "$candidate" | tr -d '\n' | shasum -a 256 | cut -c1-12)
    test "${candidate_id##*h}" = "$digest"
    test "$(jq -r .config.specialise_profile "$input")" = true
    test "$(jq -r .record.implementation.flags.specialise_profile "$candidate")" = true
done < "$out/FREEZE.tsv"
test ! -e "$out/runs"
mkdir "$out/runs"
uname -srm > "$out/runs/host.txt"
if ! sysctl -n machdep.cpu.brand_string > "$out/runs/cpu.txt" 2> "$out/runs/cpu.stderr.txt"; then
    printf 'CPU model unavailable from sysctl in this sandbox\n' > "$out/runs/cpu.txt"
fi
printf 'unisolated-macos-exploratory\n' > "$out/runs/isolation.txt"
printf 'target\trepetition\trun_id\tinput_sha256\tcandidate_sha256\tworkload_id\n' > "$out/runs/INDEX.tsv"
for target in 1 7; do
    if [ "$target" -eq 1 ]; then schedule='1 2'; else schedule='2 1'; fi
    workload_id=$(jq -r .workload_id "$archived/T$target/workload.json")
    input="$out/inputs/T${target}-f6_ic_specialise_breakdown.json"
    candidate="$out/candidates/T${target}-f6_ic_specialise_breakdown.json"
    candidate_id=$(jq -r .candidate_id "$candidate")
    for repetition in $schedule; do
        run_id="${candidate_id}W${workload_id#W}R${repetition}"
        stem="$out/runs/T${target}-R${repetition}"
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$target" "$repetition" \
            "$run_id" "$(sha "$input")" "$(sha "$candidate")" "$workload_id" \
            >> "$out/runs/INDEX.tsv"
        date -u '+%Y-%m-%dT%H:%M:%SZ' > "$stem.start-utc"
        if RAYON_NUM_THREADS=1 IC_ARTIFACT_CACHE=off IC_F2_BACKEND=cpu \
            gtimeout 180 "$binary" < "$input" > "$stem.stdout.json" \
            2> "$stem.stderr.txt"; then rc=0; else rc=$?; fi
        printf '%s\n' "$rc" > "$stem.exit-status"
        date -u '+%Y-%m-%dT%H:%M:%SZ' > "$stem.end-utc"
    done
done
shasum -a 256 "$out"/runs/*.stdout.json "$out"/runs/*.stderr.txt \
    > "$out/runs/RAW_SHA256SUMS"

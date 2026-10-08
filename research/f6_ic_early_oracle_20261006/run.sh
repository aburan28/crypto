#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
cd "$repo"
out=research/f6_ic_early_oracle_20261006
archived=research/f6_ic_geometric_closure_20261003/eight_target
binary=target/bench_bins/f6_ic_early_oracle_worker
sha() { shasum -a 256 "$1" | awk '{print $1}'; }
expected_binary=$(sed -n 's/^binary_sha256=//p' "$out/BUILD_IDENTITY.txt")
test "$(sha "$binary")" = "$expected_binary"
shasum -a 256 -c "$out/SOURCE_SHA256SUMS"
while IFS="$(printf '\t')" read -r target arm input_sha candidate_id candidate_sha; do
    [ "$target" = target ] && continue
    case "$arm" in
        early_baseline|early_candidate) input="$out/inputs/T${target}-$arm.json" ;;
        *) exit 2 ;;
    esac
    candidate="$out/candidates/$arm.json"
    test "$(sha "$input")" = "$input_sha"
    test "$(sha "$candidate")" = "$candidate_sha"
    test "$(jq -r .candidate_id "$candidate")" = "$candidate_id"
    test "$(jq -r .record.implementation.flags.active_multipliers "$candidate")" = false
    test "$(jq -r .record.implementation.flags.direct_fused_pack "$candidate")" = false
    test "$(jq -r .record.implementation.flags.support_local_stream "$candidate")" = false
    test "$(jq -r .record.implementation.flags.support_local_bitmap_columns "$candidate")" = false
    test "$(jq -r .record.implementation.flags.support_local_profile "$candidate")" = false
    test "$(jq -r .record.implementation.flags.inherited_basis_profile "$candidate")" = false
    expected_early=false
    if [ "$arm" = early_candidate ]; then expected_early=true; fi
    test "$(jq -r .record.implementation.flags.early_oracle "$candidate")" = "$expected_early"
    digest=$(jq -S -c '.record' "$candidate" | tr -d '\n' | shasum -a 256 | cut -c1-12)
    test "${candidate_id##*h}" = "$digest"
    workload="$archived/T$target/workload.json"
    test -f "$workload"
    test "$(jq -r .config.solver "$input")" = f6_ic
    test "$(jq -r .config.active_multipliers "$input")" = false
    test "$(jq -r .config.direct_fused_pack "$input")" = false
    test "$(jq -r .config.support_local_stream "$input")" = false
    test "$(jq -r .config.support_local_bitmap_columns "$input")" = false
    test "$(jq -r .config.support_local_profile "$input")" = false
    test "$(jq -r .config.inherited_basis_profile "$input")" = false
    test "$(jq -r .config.early_oracle "$input")" = "$expected_early"
done < "$out/FREEZE.tsv"
test ! -e "$out/runs"
mkdir "$out/runs"
uname -srm > "$out/runs/host.txt"
if ! sysctl -n machdep.cpu.brand_string > "$out/runs/cpu.txt" 2> "$out/runs/cpu.stderr.txt"; then
    printf 'CPU model unavailable from sysctl in this sandbox\n' > "$out/runs/cpu.txt"
fi
printf 'unisolated-macos-exploratory\n' > "$out/runs/isolation.txt"
printf 'target\trepetition\tarm\trun_id\tinput_sha256\tcandidate_sha256\tworkload_id\n' > "$out/runs/INDEX.tsv"
for target in 1 7; do
    if [ "$target" -eq 1 ]; then
        schedule='early_baseline:1 early_candidate:1 early_candidate:2 early_baseline:2'
    else
        schedule='early_candidate:1 early_baseline:1 early_baseline:2 early_candidate:2'
    fi
    workload_id=$(jq -r .workload_id "$archived/T$target/workload.json")
    for item in $schedule; do
        arm=${item%:*}
        repetition=${item#*:}
        input="$out/inputs/T${target}-$arm.json"
        candidate="$out/candidates/$arm.json"
        candidate_id=$(jq -r .candidate_id "$candidate")
        run_id="${candidate_id}W${workload_id#W}R${repetition}"
        stem="$out/runs/T${target}-R${repetition}-${arm}"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$target" "$repetition" "$arm" \
            "$run_id" "$(sha "$input")" "$(sha "$candidate")" "$workload_id" \
            >> "$out/runs/INDEX.tsv"
        date -u '+%Y-%m-%dT%H:%M:%SZ' > "$stem.start-utc"
        if RAYON_NUM_THREADS=1 IC_ARTIFACT_CACHE=off IC_F2_BACKEND=cpu \
            gtimeout 180 "$binary" < "$input" > "$stem.stdout.json" \
            2> "$stem.stderr.txt"; then
            rc=0
        else
            rc=$?
        fi
        printf '%s\n' "$rc" > "$stem.exit-status"
        date -u '+%Y-%m-%dT%H:%M:%SZ' > "$stem.end-utc"
    done
done
shasum -a 256 "$out"/runs/*.stdout.json "$out"/runs/*.stderr.txt \
    > "$out/runs/RAW_SHA256SUMS"

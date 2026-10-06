#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
cd "$repo"
out=research/f6_ic_adaptive_pair_20261005
base=research/f6_ic_geometric_closure_20261003/eight_target
binary=target/bench_bins/f6_ic_adaptive_pair_worker
sha() { shasum -a 256 "$1" | awk '{print $1}'; }
expected_binary=$(sed -n 's/^binary_sha256=//p' "$out/BUILD_IDENTITY.txt")
test "$(sha "$binary")" = "$expected_binary"
shasum -a 256 -c "$out/SOURCE_SHA256SUMS"
while IFS="$(printf '\t')" read -r target arm input_sha candidate_id candidate_sha; do
    [ "$target" = target ] && continue
    case "$arm" in
        f6_ic_pair) input="$out/inputs/T${target}-$arm.json" ;;
        inherited_f4|f6_ic) input="$base/T$target/input-$arm.json" ;;
        *) exit 2 ;;
    esac
    candidate="$out/candidates/$arm.json"
    test "$(sha "$input")" = "$input_sha"
    test "$(sha "$candidate")" = "$candidate_sha"
    test "$(jq -r .candidate_id "$candidate")" = "$candidate_id"
    digest=$(jq -S -c '.record' "$candidate" | tr -d '\n' | shasum -a 256 | cut -c1-12)
    test "${candidate_id##*h}" = "$digest"
    workload="$base/T$target/workload.json"
    test -f "$workload"
    test "$(jq -r .config.solver "$input")" = "$arm"
done < "$out/FREEZE.tsv"
test ! -e "$out/runs"
mkdir "$out/runs"
uname -srm > "$out/runs/host.txt"
sysctl -n machdep.cpu.brand_string > "$out/runs/cpu.txt"
printf 'unisolated-macos-exploratory\n' > "$out/runs/isolation.txt"
printf 'target\trepetition\tarm\trun_id\tinput_sha256\tcandidate_sha256\tworkload_id\n' > "$out/runs/INDEX.tsv"
for target in 1 7; do
    if [ "$target" -eq 1 ]; then
        schedule='inherited_f4:1 f6_ic:1 f6_ic_pair:1 f6_ic_pair:2 f6_ic:2 inherited_f4:2'
    else
        schedule='f6_ic_pair:1 f6_ic:1 inherited_f4:1 inherited_f4:2 f6_ic:2 f6_ic_pair:2'
    fi
    workload_id=$(jq -r .workload_id "$base/T$target/workload.json")
    for item in $schedule; do
        arm=${item%:*}
        repetition=${item#*:}
        case "$arm" in
            f6_ic_pair) input="$out/inputs/T${target}-$arm.json" ;;
            *) input="$base/T$target/input-$arm.json" ;;
        esac
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

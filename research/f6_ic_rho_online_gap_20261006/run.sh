#!/bin/sh
set -eu
repo=$(CDPATH= cd -- "$(dirname -- "$0")/../.." && pwd)
cd "$repo"
out=research/f6_ic_rho_online_gap_20261006
binary=target/bench_bins/f6_ic_compact_refuted_worker
sha() { shasum -a 256 "$1" | awk '{print $1}'; }
test "$(sha "$binary")" = "$(sed -n 's/^binary_sha256=//p' "$out/BUILD_IDENTITY.txt")"
shasum -a 256 -c "$out/SOURCE_SHA256SUMS"
shasum -a 256 -c "$out/FROZEN_SHA256SUMS"
test ! -e "$out/runs"
while IFS="$(printf '\t')" read -r target workload candidate ic_sha rho_sha point_sha; do
    [ "$target" = target ] && continue
    ic="$out/inputs/T$target-ic.json"
    rho="$out/inputs/T$target-rho.json"
    test "$(sha "$ic")" = "$ic_sha"
    test "$(sha "$rho")" = "$rho_sha"
    test "$(jq -r .candidate_id "$out/candidates/f6_ic.json")" = "$candidate"
    test "$(jq -r .workload_id "$out/workloads/T$target.json")" = "$workload"
    test "$(jq -S -c '.public_targets[0]' "$rho" | shasum -a 256 | cut -d ' ' -f 1)" = "$point_sha"
    test "$(jq -r .mode "$ic")" = ic
    test "$(jq -r .mode "$rho")" = rho
    test "$(jq -r .config.rho_parallel_walks "$rho")" = 1
    test "$(jq -r .config.max_trials "$rho")" = 65536
done < "$out/FREEZE.tsv"
mkdir "$out/runs"
uname -srm > "$out/runs/host.txt"
if ! sysctl -n machdep.cpu.brand_string > "$out/runs/cpu.txt" 2> "$out/runs/cpu.stderr.txt"; then
    printf 'CPU model unavailable from sysctl in this sandbox\n' > "$out/runs/cpu.txt"
fi
printf 'unisolated-macos-exploratory\n' > "$out/runs/isolation.txt"
printf 'target\trepetition\tarm\trun_id\tworkload_id\tinput_sha256\tpoint_sha256\n' > "$out/runs/INDEX.tsv"
for target in 1 7; do
    workload=$(jq -r .workload_id "$out/workloads/T$target.json")
    point_sha=$(jq -S -c '.public_targets[0]' "$out/inputs/T$target-ic.json" | shasum -a 256 | cut -d ' ' -f 1)
    candidate=$(jq -r .candidate_id "$out/candidates/f6_ic.json")
    if [ "$target" -eq 1 ]; then
        schedule='ic:1 rho:1 rho:2 ic:2'
    else
        schedule='rho:1 ic:1 ic:2 rho:2'
    fi
    for item in $schedule; do
        arm=${item%:*}
        repetition=${item#*:}
        input="$out/inputs/T$target-$arm.json"
        if [ "$arm" = ic ]; then
            run_id="${candidate}W${workload#W}R${repetition}"
        else
            run_id="rho-${workload}-R${repetition}"
        fi
        stem="$out/runs/T${target}-R${repetition}-${arm}"
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
          "$target" "$repetition" "$arm" "$run_id" "$workload" "$(sha "$input")" "$point_sha" \
          >> "$out/runs/INDEX.tsv"
        date -u '+%Y-%m-%dT%H:%M:%SZ' > "$stem.start-utc"
        if RAYON_NUM_THREADS=1 IC_ARTIFACT_CACHE=off IC_F2_BACKEND=cpu \
            gtimeout 180 "$binary" < "$input" > "$stem.stdout.json" 2> "$stem.stderr.txt"; then
            rc=0
        else
            rc=$?
        fi
        printf '%s\n' "$rc" > "$stem.exit-status"
        date -u '+%Y-%m-%dT%H:%M:%SZ' > "$stem.end-utc"
    done
done
shasum -a 256 "$out"/runs/*.stdout.json "$out"/runs/*.stderr.txt > "$out/runs/RAW_SHA256SUMS"

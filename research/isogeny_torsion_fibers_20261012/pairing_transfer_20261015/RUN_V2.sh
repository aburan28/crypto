#!/bin/sh
# Thin, fixed-command orchestration only. The native Rust example performs
# all arithmetic, selection, solving, and verification.
set -eu

study=research/isogeny_torsion_fibers_20261012/pairing_transfer_20261015
bin=target/release/examples/torsion_pairing_transfer
input=$study/INPUTS.jsonl
out=$study/runs_v2

check_sha() {
    actual=$(shasum -a 256 "$1" | awk '{print $1}')
    if [ "$actual" != "$2" ]; then
        printf 'SHA mismatch: %s\n' "$1" >&2
        exit 1
    fi
}

check_sha examples/torsion_pairing_transfer.rs 382ea60f8746c5145ebffe6d25f6809a520b09974f9f5e294da9d96b96d82de7
check_sha "$input" ef5cd2ac1d8a2cb9aa3b7edc6f5b92408fafc00e7253b97dddc0d3555998ddc5
check_sha "$study/SELECTION_61.json" 911d9a45e554812a26e554965f117fa4397bc6e2bb9317b2dc8ba9c65d28c9c4
check_sha "$study/SELECTION_271.json" 958e910e8ad22092ffbef9a97b9a4636382a261273db0bd4e971b243dab210bc
check_sha "$bin" 9daceafea0dce4811299ff87a6f4ed3818dc1b49ed5410dffb2133fa56af804d

mkdir -p "$out"
{
    date -u '+start_utc=%Y-%m-%dT%H:%M:%SZ'
    git rev-parse HEAD | sed 's/^/commit=/'
    uname -a | sed 's/^/host=/'
    rustc --version | sed 's/^/rustc=/'
    printf 'isolation=L0_shared_host_sequential\n'
    printf 'source_sha256=382ea60f8746c5145ebffe6d25f6809a520b09974f9f5e294da9d96b96d82de7\n'
    printf 'binary_sha256=9daceafea0dce4811299ff87a6f4ed3818dc1b49ed5410dffb2133fa56af804d\n'
} > "$out/RUN_META.txt"
printf 'name\tcommand\texit\n' > "$out/RUN_STATUS.tsv"

run_one() {
    name=$1
    shift
    code=0
    "$@" > "$out/$name.json" 2> "$out/$name.stderr" || code=$?
    printf '%s\t%s\t%d\n' "$name" "$*" "$code" >> "$out/RUN_STATUS.tsv"
    if [ "$code" -ne 0 ]; then
        printf 'Stopped at %s (exit %d); raw output retained.\n' "$name" "$code" >&2
        exit "$code"
    fi
}

for cell in 0 1; do
    if [ "$cell" -eq 0 ]; then selection=$study/SELECTION_61.json; else selection=$study/SELECTION_271.json; fi
    for round in 0 1 2; do
        prefix="c${cell}_r${round}"
        run_one "${prefix}_pairing" "$bin" run "$input" "$cell" "$selection"
        run_one "${prefix}_verify" "$bin" verify "$input" "$cell" "$out/${prefix}_pairing.json"
        run_one "${prefix}_bilinear" "$bin" bilinear "$input" "$cell" "$selection"
        run_one "${prefix}_rho" "$bin" rho "$input" "$cell" "$round"
        run_one "${prefix}_rho_aa" "$bin" rho "$input" "$cell" "$round"
    done
done

date -u '+end_utc=%Y-%m-%dT%H:%M:%SZ' >> "$out/RUN_META.txt"

#!/usr/bin/env bash
# Model-check StrongRho.tla with TLC under every *.cfg in this directory.
#
#   ./run_tlc.sh            # all configs
#   ./run_tlc.sh a.cfg ...  # only these
#
# Each cfg names its expected verdict on a line "\* EXPECT: ...":
#   pass                 TLC must finish with "No error has been found"
#   violate <Invariant>  TLC must report "Invariant <Invariant> is violated"
#   violate temporal     TLC must report a temporal-property violation
# Output of each run goes to out/<config>.txt and a table to out/SUMMARY.txt.
# Exits nonzero when any config does not meet its expectation.
set -u -o pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CACHE="${TLA_CACHE:-/tmp/claude-0/tla}"
JAR="$CACHE/tla2tools.jar"
JAR_URL="https://github.com/tlaplus/tlaplus/releases/download/v1.8.0/tla2tools.jar"
# Pinned after the first download on 2026-10-16 (the jar reports
# "TLC2 Version 2026.10.15.145915 (rev: 932d35e)").
JAR_SHA256="15d103b7e5617780bcc904c0af12f1ad242e2b345c0301136c5be4588f17edac"
WORKERS="${TLC_WORKERS:-2}"   # shared machine: at most 2
# TLC spills its state queue and fingerprints to -metadir; keep that in RAM
# (tmpfs) when there is one, since the shared disk can fill up mid-run.
if [[ -d /dev/shm && -w /dev/shm ]]; then META_ROOT="${TLC_META:-/dev/shm}"; else META_ROOT="${TLC_META:-$CACHE}"; fi

mkdir -p "$CACHE" "$HERE/out"
if [[ ! -f "$JAR" ]]; then
    echo "downloading $JAR_URL"
    curl -fsSL -o "$JAR.part" "$JAR_URL" && mv "$JAR.part" "$JAR" || { echo "download failed" >&2; exit 2; }
fi
got="$(sha256sum "$JAR" | cut -d' ' -f1)"
if [[ "$got" != "$JAR_SHA256" ]]; then
    echo "SHA-256 mismatch for $JAR: got $got, want $JAR_SHA256" >&2
    exit 2
fi

cd "$HERE"
if [[ $# -gt 0 ]]; then cfgs=("$@"); else cfgs=( $(ls *.cfg | sort) ); fi

summary="$HERE/out/SUMMARY.txt"
[[ $# -gt 0 ]] || : > "$summary"
fails=0
for cfg in "${cfgs[@]}"; do
    name="$(basename "$cfg" .cfg)"
    expect="$(sed -n 's/^\\\* EXPECT: *//p' "$cfg" | head -1)"
    [[ -n "$expect" ]] || { echo "$cfg: no EXPECT line" >&2; fails=$((fails+1)); continue; }
    out="$HERE/out/$name.txt"
    meta="$(mktemp -d "$META_ROOT/strongrho-tlc.XXXXXX")"
    start=$(date +%s)
    java -XX:+UseParallelGC -cp "$JAR" tlc2.TLC -workers "$WORKERS" -metadir "$meta" \
        -config "$cfg" StrongRho.tla 2>&1 | grep -v '^Picked up JAVA_TOOL_OPTIONS' > "$out"
    secs=$(( $(date +%s) - start ))
    rm -rf "$meta"

    distinct="$(sed -n 's/.* \([0-9,]*\) distinct states found.*/\1/p' "$out" | tail -1)"
    depth="$(sed -n 's/The depth of the complete state graph search is \([0-9]*\)\./\1/p' "$out" | tail -1)"
    if grep -q '^Model checking completed. No error has been found.' "$out"; then
        verdict="pass"
    elif v="$(sed -n 's/^Error: Invariant \([A-Za-z0-9_]*\) is violated.*/\1/p' "$out" | head -1)"; [[ -n "$v" ]]; then
        verdict="violate $v"
    elif grep -q 'Temporal properties were violated' "$out"; then
        verdict="violate temporal"
    elif grep -q 'Error: Deadlock reached' "$out"; then
        verdict="deadlock"
    else
        verdict="error"
    fi
    if [[ "$verdict" == "$expect" ]]; then ok="OK"; else ok="FAIL"; fails=$((fails+1)); fi
    line="$(printf '%-34s expect=%-26s got=%-26s distinct=%-11s depth=%-4s %4ss %s' \
        "$name" "$expect" "$verdict" "${distinct:--}" "${depth:--}" "$secs" "$ok")"
    echo "$line"
    echo "$line" >> "$summary"
done
echo "failures: $fails"
[[ $fails -eq 0 ]]

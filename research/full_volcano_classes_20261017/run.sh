#!/usr/bin/env bash
# Thin native-CLI orchestration. Each case owns a separate SQLite database.
set -euo pipefail
bin=${1:?isogeny-class executable}
out=${2:?output directory}
name=${3:?two, three, or mixed}
case "$name" in two|three|mixed) ;; *) exit 2 ;; esac
study=$(cd -- "$(dirname -- "$0")" && pwd)
mkdir -p -- "$out"
code=0
timeout 7200 "$bin" build --curve "$study/inputs/$name.json" --out "$out/$name.sqlite" --steps all > "$study/evidence/$name-build.json" 2> "$study/evidence/$name-build.stderr" || code=$?
printf '%s\n' "$code" > "$study/evidence/$name-build.exit"
[ "$code" = 0 ] || exit "$code"
sha256sum "$out/$name.sqlite" > "$study/evidence/$name-before-replay.sha256"
cp -- "$out/$name.sqlite" "$out/$name-replay.sqlite"
code=0
timeout 7200 "$bin" verify --input "$out/$name-replay.sqlite" > "$study/evidence/$name-verify.json" 2> "$study/evidence/$name-verify.stderr" || code=$?
printf '%s\n' "$code" > "$study/evidence/$name-verify.exit"
[ "$code" = 0 ] || exit "$code"
sha256sum -c "$study/evidence/$name-before-replay.sha256"

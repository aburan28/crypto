#!/usr/bin/env bash
# E3 (Amendment 2): random base, 8 draws, rand and sat. Thin orchestration.
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"; root="$(cd "$here/../../.." && pwd)"
cd "$root"; cargo build --release --example hky_last_fall_degree
bin="$root/target/release/examples/hky_last_fall_degree"
git rev-parse HEAD > "$here/commit"; sha256sum "$bin" | cut -d' ' -f1 > "$here/sha256"
cell() { "$bin" cell --n "$1" --family random --target "$2" --draws 8 --seed 20261012 --cmax 6 > "$here/n$1_$2.jsonl" 2> "$here/n$1_$2.err" && echo "$(date -u +%FT%TZ) ok n=$1 $2" >> "$here/cells.log" || echo "$(date -u +%FT%TZ) FAIL n=$1 $2" >> "$here/cells.log"; }
export -f cell; export bin here
# small cells four at a time, then n = 32 two at a time
printf '%s\n' "18 rand" "18 sat" "20 rand" "20 sat" "22 rand" "22 sat" "24 rand" "24 sat" "26 rand" "26 sat" "28 rand" "28 sat" | xargs -P 4 -L 1 bash -c 'cell $0 $1'
printf '%s\n' "32 rand" "32 sat" | xargs -P 2 -L 1 bash -c 'cell $0 $1'

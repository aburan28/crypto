#!/bin/zsh
set -euo pipefail

study=research/isogeny_torsion_fibers_20261012/different_j_order_screen_20261015
[[ -f Cargo.toml && -f "$study/RUN_SPEC_V1.json" ]] || {
  print -u2 'run from the crypto repository root after committing RUN_SPEC_V1.json'
  exit 2
}
for output in RAW_V1.jsonl RAW_V1.stderr RUN_V1.exit RUN_V1.started_utc RUN_V1.finished_utc REPLAY_V1.json REPLAY_V1.stderr REPLAY_V1.exit; do
  [[ ! -e "$study/$output" ]] || {
    print -u2 "refusing to overwrite $study/$output"
    exit 2
  }
done

file_sha() { shasum -a 256 "$1" | cut -d ' ' -f 1 }
[[ "$(file_sha examples/different_j_order_screen.rs)" == ed0842099f99f902dbb90dcc8b0c512db090c8ecbf31fd9ca873c43aee3d8df0 ]]
[[ "$(file_sha examples/different_j_order_check.rs)" == 402d1ed6326a5f3c69711231d702e06271495e7849bbc84e35bd4094196cbbb4 ]]
[[ "$(file_sha target/release/examples/different_j_order_screen)" == 5395ebbb7ca1b3b2befdeeb016fd99e0ab19287ee36650cdbab92beef7491536 ]]
[[ "$(file_sha target/release/examples/different_j_order_check)" == 5d8b3aa2b5b8c9b4e0a4fb7b22446d61167721a00f8862ee055b57a711ec5c3e ]]

date -u '+%Y-%m-%dT%H:%M:%SZ' > "$study/RUN_V1.started_utc"
set +e
/opt/homebrew/bin/gtimeout --signal=TERM --kill-after=5s 120s target/release/examples/different_j_order_screen screen > "$study/RAW_V1.jsonl" 2> "$study/RAW_V1.stderr"
run_exit=$?
set -e
print -r -- "$run_exit" > "$study/RUN_V1.exit"
date -u '+%Y-%m-%dT%H:%M:%SZ' > "$study/RUN_V1.finished_utc"
[[ "$run_exit" == 0 ]] || exit "$run_exit"

set +e
target/release/examples/different_j_order_check "$study/RAW_V1.jsonl" > "$study/REPLAY_V1.json" 2> "$study/REPLAY_V1.stderr"
replay_exit=$?
set -e
print -r -- "$replay_exit" > "$study/REPLAY_V1.exit"
exit "$replay_exit"

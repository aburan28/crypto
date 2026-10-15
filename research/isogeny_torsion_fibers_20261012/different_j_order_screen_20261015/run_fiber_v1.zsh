#!/bin/zsh
set -euo pipefail

study=research/isogeny_torsion_fibers_20261012/different_j_order_screen_20261015
[[ -f Cargo.toml && -f "$study/FIBER_SPEC_V1.json" ]] || exit 2
for output in FIBER_V1.json FIBER_V1.stderr FIBER_V1.exit FIBER_REPLAY_V1.json FIBER_REPLAY_V1.stderr FIBER_REPLAY_V1.exit; do
  [[ ! -e "$study/$output" ]] || {
    print -u2 "refusing to overwrite $study/$output"
    exit 2
  }
done
file_sha() { shasum -a 256 "$1" | cut -d ' ' -f 1 }
[[ "$(file_sha "$study/FIBER_PROTOCOL_V1.md")" == 535ee294315eba03353495944faa378b949c64a5e8d92221160f62f82b2240c5 ]]
[[ "$(file_sha "$study/RAW_V1.jsonl")" == 31238914df3663805efb88c3c89c420b4ca42b875668da0cfc954cd1084171b4 ]]
[[ "$(file_sha "$study/CENSUS_V1.jsonl")" == 50f2b79d9811c72ccdacc316edf53b657b0d3a7911a332b64cc70749d315b925 ]]
[[ "$(file_sha examples/different_j_order_fiber.rs)" == 8f2dd30398140c8f5fc746fb0e059b0527ff175f64aa5fa85e8c26dcd316e0b2 ]]
[[ "$(file_sha examples/different_j_order_fiber_check.rs)" == 8a20e2afd09ac8172d461eef63594e2ee2be95870cc64b9d39cd5fc9af2bc021 ]]
[[ "$(file_sha target/release/examples/different_j_order_fiber)" == 5582b21f40126be75a3cd19f01029c366d7e57a1fb1dc6871b61e4ed505f1004 ]]
[[ "$(file_sha target/release/examples/different_j_order_fiber_check)" == a4a9ae7bbff128d7f7a1f810cffaf02bf229e7f3b0374bca508f31d8384a440a ]]

set +e
/opt/homebrew/bin/gtimeout --signal=TERM --kill-after=5s 120s target/release/examples/different_j_order_fiber audit > "$study/FIBER_V1.json" 2> "$study/FIBER_V1.stderr"
fiber_exit=$?
set -e
print -r -- "$fiber_exit" > "$study/FIBER_V1.exit"
[[ "$fiber_exit" == 0 ]] || exit "$fiber_exit"

set +e
/opt/homebrew/bin/gtimeout --signal=TERM --kill-after=5s 120s target/release/examples/different_j_order_fiber_check "$study/FIBER_V1.json" > "$study/FIBER_REPLAY_V1.json" 2> "$study/FIBER_REPLAY_V1.stderr"
replay_exit=$?
set -e
print -r -- "$replay_exit" > "$study/FIBER_REPLAY_V1.exit"
exit "$replay_exit"

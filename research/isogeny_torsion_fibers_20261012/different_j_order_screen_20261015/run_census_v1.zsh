#!/bin/zsh
set -euo pipefail

study=research/isogeny_torsion_fibers_20261012/different_j_order_screen_20261015
[[ -f Cargo.toml && -f "$study/CENSUS_SPEC_V1.json" ]] || exit 2
for output in CENSUS_V1.jsonl CENSUS_V1.stderr CENSUS_V1.exit; do
  [[ ! -e "$study/$output" ]] || {
    print -u2 "refusing to overwrite $study/$output"
    exit 2
  }
done
file_sha() { shasum -a 256 "$1" | cut -d ' ' -f 1 }
[[ "$(file_sha examples/different_j_order_census.rs)" == 6086a7b8d27f94a870d70749c5478357fe37f23c2826767e36e06b53fc6ad61c ]]
[[ "$(file_sha target/release/examples/different_j_order_census)" == 1112577f003d2326f607a62a8142d55623172cd5229b2e8b0018876f8b70b974 ]]
[[ "$(file_sha "$study/RAW_V1.jsonl")" == 31238914df3663805efb88c3c89c420b4ca42b875668da0cfc954cd1084171b4 ]]

set +e
/opt/homebrew/bin/gtimeout --signal=TERM --kill-after=5s 120s target/release/examples/different_j_order_census "$study/RAW_V1.jsonl" > "$study/CENSUS_V1.jsonl" 2> "$study/CENSUS_V1.stderr"
census_exit=$?
set -e
print -r -- "$census_exit" > "$study/CENSUS_V1.exit"
exit "$census_exit"

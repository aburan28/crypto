#!/bin/zsh
set -euo pipefail

study_dir=research/isogeny_torsion_fibers_20261012/different_j_full_image_20261015
public_input="$study_dir/PUBLIC_V2.json"
binary=target/release/examples/nonisomorphic_extension_isogeny
record_dir="$study_dir/raw-v2"
combined="$study_dir/RUN_V2.jsonl"

if [[ -e "$combined" || -e "$record_dir" ]]; then
  print -u2 'V2 output already exists; preserve it and use a new version'
  exit 3
fi
mkdir -p "$record_dir"
: > "$combined"

for round in 0 1 2; do
  case "$round" in
    0) arms=(source isogeny-image destination-direct full-image rho-negation rho-negation-aa) ;;
    1) arms=(isogeny-image destination-direct full-image rho-negation rho-negation-aa source) ;;
    2) arms=(destination-direct full-image rho-negation rho-negation-aa source isogeny-image) ;;
  esac
  for arm in "${arms[@]}"; do
    stem="$record_dir/round-${round}-${arm}"
    if [[ "$arm" == rho-negation ]]; then
      if "$binary" rho "$public_input" "$round" reference > "$stem.stdout" 2> "$stem.stderr"; then
        exit_code=0
      else
        exit_code=$?
      fi
    elif [[ "$arm" == rho-negation-aa ]]; then
      if "$binary" rho "$public_input" "$round" aa > "$stem.stdout" 2> "$stem.stderr"; then
        exit_code=0
      else
        exit_code=$?
      fi
    else
      if "$binary" run "$public_input" "$round" "$arm" > "$stem.stdout" 2> "$stem.stderr"; then
        exit_code=0
      else
        exit_code=$?
      fi
    fi
    print -- "$exit_code" > "$stem.exit"
    if (( exit_code != 0 )); then
      print -u2 "run failed: round=$round arm=$arm exit=$exit_code"
      exit "$exit_code"
    fi
    cat "$stem.stdout" >> "$combined"
  done
done

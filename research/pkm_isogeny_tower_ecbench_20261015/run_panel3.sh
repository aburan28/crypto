#!/bin/bash
# Panel 3 of PROTOCOL.md: the auxiliary-curve search's cost distribution.
# `ecbench fb` builds prime-isogeny:t=T,seed=S on the 24-bit curve of each
# panel, t = 4 … 14, S = 0 … 15, and this keeps one line per build:
# curve, t, seed, fb_id, abscissae, auxiliary curves tried, charged
# additions and doublings (the whole build: search + tower).  The dumps
# themselves go to a scratch directory and are deleted.
# From the repository root:
#   ECBENCH=path/to/ecbench research/pkm_isogeny_tower_ecbench_20261015/run_panel3.sh
set -u
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
bin=${ECBENCH:-$root/target/release/ecbench}
tmp=$(mktemp -d)
out=$here/tables/panel3.tsv
mkdir -p "$here/tables"
printf 'curve\tt\tseed\tfb_id\tabscissae\taux_curves_tried\tbuild_adds\tbuild_doubles\n' > "$out"
for name in p24 q24; do
  case $name in p24) line=$(sed -n 3p "$here/curves-panel1.jsonl");; q24) line=$(sed -n 3p "$here/curves-panel2.jsonl");; esac
  for t in 4 5 6 7 8 9 10 11 12 13 14; do
    for s in $(seq 0 15); do
      nice -n 10 prlimit --as=2684354560 -- "$bin" fb --curve "$line" --factor-base "prime-isogeny:t=$t,seed=$s" --out "$tmp/fb.json" > /dev/null 2>&1 || { printf '%s\t%s\t%s\tFAILED\n' $name $t $s >> "$out"; continue; }
      jq -r --arg n $name --arg t $t --arg s $s '[$n, $t, $s, .factor_base.fb_id, .factor_base.abscissae, (.factor_base.description | capture("(?<c>[0-9]+) auxiliary curves tried").c), .build_adds, .build_doubles] | @tsv' "$tmp/fb.json" >> "$out"
      rm -f "$tmp/fb.json"
    done
  done
done
rmdir "$tmp"

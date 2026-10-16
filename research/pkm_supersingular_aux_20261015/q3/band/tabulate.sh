#!/bin/bash
# Lay out the ecbench band reports (band/*.band.json, ecbench.band/v1): each
# arm's mean priced S at w = 0, 0.1, 1, and for the study's pairs the ratio of
# those means (B / A) with ecbench band's own verdict (division only).
set -eu
here=$(cd "$(dirname "$0")" && pwd)
cd "$here"
{
  printf 'curve\tarm\tS_w0\tS_w0.1\tS_w1\n'
  for n in s20 s22 s24; do
    jq -r --arg n "$n" '.rows[] | [$n, .arm, .mean_s[0], .mean_s[1], .mean_s[2]] | @tsv' $n.band.json
  done
} > band-arms.tsv
{
  printf 'curve\tB\tA\tratio_w0\tratio_w0.1\tratio_w1\tband_verdict\n'
  for n in s20 s22 s24; do
    for pair in ic-ss-rootfind:ic-plain-ss ic-ord-rootfind:ic-plain-ord ic-ss-rootfind:ic-ord-rootfind ic-plain-ss:ic-plain-ord ic-ss-rootfind:ic-ss-mitm ic-ss-rootfind:rho-neg ic-ord-rootfind:rho-neg rho-neg-aa:rho-neg; do
      b=${pair%%:*}; a=${pair##*:}
      jq -r --arg n "$n" --arg a "$a" --arg b "$b" '
        (.rows | map({(.arm): .mean_s}) | add) as $m |
        (.ranks[] | select((.a == $a and .b == $b) or (.a == $b and .b == $a)) | .verdict) as $v |
        [$n, $b, $a, ($m[$b][0] / $m[$a][0]), ($m[$b][1] / $m[$a][1]), ($m[$b][2] / $m[$a][2]), $v] | @tsv' $n.band.json
    done
  done
} > band-ratios.tsv

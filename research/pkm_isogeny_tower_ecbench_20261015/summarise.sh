#!/bin/bash
# Lays out recorded values (no new measurement): per-session ratios from the
# saved `ecbench compare` files, per-arm mean S from the sealed records, and
# the isogeny arms' factor-base phase (the charged auxiliary-curve search)
# against their total, into tables/summary-*.tsv.  jq only.
set -eu
here=$(cd "$(dirname "$0")" && pwd)
cd "$here"
out=tables/summary-ratios.tsv
printf 'session\tA\tB\tratio_B_over_A\tlo95\thi95\tpairs\tstatus\tbounded\n' > $out
for d in sessions/*/; do
  n=$(basename $d)
  for f in $d/comparisons/*.json; do
    jq -r --arg n $n '[$n, .a.arm, .b.arm, .ops.ratio_b_over_a, .ops.ci95[0], .ops.ci95[1], .ops.pairs, .ops.status, .ops.bounded] | @tsv' $f >> $out
  done
done
out=tables/summary-arms.tsv
printf 'session\tarm\truns\tverified\tmean_S\tmean_total_gae\tsqrt_r\tfb_gae\tmean_fb_share\taux_curves_tried\taux_search_gae\taux_count_gae\taux_build_gae\tabscissae\tfb_id\n' > $out
for d in sessions/*/; do
  n=$(basename $d)
  jq -s -r --arg n $n '
    map(select(.warmup == false)) | group_by(.arm)[] |
    . as $runs | ($runs[0]) as $r0 |
    ([$runs[] | .phases[] | select(.name == "factor_base") | .gae] ) as $fb |
    ($r0.phases // [] | map(select(.name == "factor_base"))[0].native // {}) as $nat |
    [$n, $r0.arm, ($runs|length), ($runs | map(select(.outcome.status == "verified")) | length),
     ($runs | map(.cost.s) | add / length),
     ($runs | map(.cost.total_gae) | add / length),
     ($r0.workload.curve.r | sqrt),
     (if ($fb|length) > 0 then ($fb | add / length) else 0 end),
     (if ($fb|length) > 0 then ([$runs[] | ((.phases[] | select(.name == "factor_base") | .gae) / .cost.total_gae)] | add / length) else 0 end),
     ($nat.aux_curves_tried // ""), ($nat.aux_search_gae // ""), ($nat.aux_count_gae // ""), ($nat.aux_tower_build_gae // ""),
     ($r0.factor_base.abscissae // ""), ($r0.factor_base.fb_id // "")] | @tsv' $d/records.jsonl >> $out
done
# Panel 3: per (curve, t) mean / min / max of auxiliary curves tried over the
# sixteen seeds, the mean charged GAE of the whole build and per curve tried,
# and the pooled C = Σ tried / Σ 2^t over t = 8 … 14 (both curves).
out=tables/panel3-summary.tsv
printf 'curve\tt\tseeds\tmean_tried\tmean_tried_over_2t\tmin_tried\tmax_tried\tmean_charged_gae\tgae_per_curve_tried\n' > $out
jq -R -s -r 'split("\n")[1:] | map(select(length>0) | split("\t")) | group_by([.[0], (.[1]|tonumber)]) | map({c: .[0][0], t: (.[0][1]|tonumber), n: length, tried: (map(.[5]|tonumber)), gae: (map((.[6]|tonumber)+(.[7]|tonumber)))}) | sort_by(.c, .t)[] | [.c, .t, .n, (.tried|add/length), ((.tried|add/length) / pow(2; .t)), (.tried|min), (.tried|max), (.gae|add/length), ((.gae|add)/(.tried|add))] | @tsv' tables/panel3.tsv >> $out
jq -R -s -r 'split("\n")[1:] | map(select(length>0) | split("\t") | select((.[1]|tonumber) >= 8)) | "pooled C (t = 8 … 14, both curves, \(length) builds) = \((map(.[5]|tonumber)|add) / (map(pow(2; .[1]|tonumber))|add))"' tables/panel3.tsv > tables/panel3-pooled.txt

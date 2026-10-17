#!/bin/sh
. "$(dirname "$0")/env.sh"
cd $STUDY/calibration
# rule: largest b (of 18..26) at which the best verified IC config (lowest mean S) has mean wall <= 0.5 s and every IC run verifies
for B in 18 20 22 24 26; do
  jq -c --arg b $B '[.[]|select(.arm|startswith("ic-"))|select(.verify_rate==1)] | sort_by(.mean_S) | .[0] | {b:($b|tonumber),best:.arm,mean_S,mean_wall_s,cv_S_per_solve:.cv_S_over_records,ok:(.mean_wall_s<=0.5)}' sessions/b$B.stats.json
done > rule_table.jsonl
BSEL=$(jq -s 'map(select(.ok))|map(.b)|max' rule_table.jsonl)
BEST=$(jq -s -r --argjson b $BSEL 'map(select(.b==$b))[0].best' rule_table.jsonl)
jq -n --argjson b $BSEL --arg best "$BEST" --slurpfile st sessions/b$BSEL.stats.json --slurpfile cfg configs.json \
 --arg ca_commit "$(git -C /home/user/cryptanalysis rev-parse HEAD)" --arg ca_base "$(git -C /home/user/cryptanalysis rev-parse c773354)" \
 --arg crypto_commit "$(git -C /home/user/crypto rev-parse HEAD)" --arg vr_commit "$(git -C $SP/vr rev-parse HEAD)" \
 --arg ecb_sha "$($ECB host | awk '/^binary/{print $2}')" \
 --argjson rule "$(jq -s . rule_table.jsonl)" '
 ($cfg[0].configs|map(select(.id==$best))[0]) as $c |
 ($st[0]|map(select(.arm==$best))[0]) as $ic |
 def n($cv;$h): (pow(1.96*$cv/$h;2)|ceil);
 def nn($cv): {"h10":n($cv;0.10),"h5":n($cv;0.05),"h3":n($cv;0.03)};
 ($st[0]|map(select(.arm=="rho-neg"))[0]) as $rho | ($st[0]|map(select(.arm=="bsgs"))[0]) as $bs | ($st[0]|map(select(.arm=="bsgs-neg"))[0]) as $bn |
 {schema:"volcano-census.frozen/v1", frozen_after:"calibration sessions b=18..26, before any census search or Phase B run",
  b:$b, ic:{id:$c.id,factor_base:$c.factor_base,oracle:$c.oracle,selection:"lowest mean S among verified IC configs at the largest b whose best IC arm has mean wall <= 0.5 s; rule_table lists every b",rule_table:$rule},
  cv_estimates_at_b:{ic_per_solve:$ic.cv_S_over_records, rho_neg_per_solve:$rho.cv_S_over_records, bsgs_textbook_per_solve:$bs.cv_S_over_records, bsgs_negation_per_solve:$bn.cv_S_over_records,
    note:"sd/mean of S over the 24 measured records (8 workloads x 3 rounds) of the calibration root at b; PLAN 5 placeholders were rho 0.52, bsgs 0.19, ic 1.0. One root, 24 solves each: wide intervals; ICC not estimated (one curve)."},
  n_per_level_for_relative_half_width:{ic:nn($ic.cv_S_over_records), rho_neg:nn($rho.cv_S_over_records), bsgs_textbook:nn($bs.cv_S_over_records), formula:"n = (1.96 CV / h)^2, ICC = 0 assumed (m solves per curve; true n_eff smaller)",
     rho_floor_rule:"PLAN 5: at least 2000 rho solves per level (kept; it exceeds the table)"},
  seeds:{calibration:"calibration/seeds.json", search:"census/search_plan.json", census:"census/seeds.json"},
  tools:{cryptanalysis_base:$ca_base, cryptanalysis_branch_head:$ca_commit, cryptanalysis_branch:"claude/volcano-census-export-20261017", crypto_branch_head_at_freeze:$crypto_commit, volcano_rho_commit:$vr_commit, volcano_rho_patch:"volcano-rho-min-n.patch", ecbench_binary_sha256:$ecb_sha}}' > ../FROZEN.json
cat ../FROZEN.json

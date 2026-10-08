#!/bin/bash
# stages.sh RUNS: each scan stage in ns a scanned summand (the stage's
# counter cycles over the counter's rate, over the process's summands),
# and the admitted keys a summand, per curve: the mean over its clean
# processes. A process its isolation record marks contended is left out.
# The run tree keeps the suite's directory names; the table names curves
# by their ICV1 slugs.
set -eu
RUNS=$1
slug() {
  case $1 in
    k0n53) echo icv1-f2m53-tm56619371-dac20a85 ;;
    k1n59) echo icv1-f2m59-tm943548413-98844ecc ;;
    k0n61) echo icv1-f2m61-t158598901-ab42b6c5 ;;
    *) echo "unknown size $1" >&2; exit 1 ;;
  esac
}
printf "curve\tprocesses\tsubtract_ns\tkey_ns\tfilter_ns\tadmitted_ns\tscan_ns\ttrial_ns\tadmitted_a_summand\n"
for f in "$RUNS"/*/*/r*.price.json; do
  size=$(echo "$f" | sed 's#.*/runs/##; s#/.*##')
  iso=${f%.price.json}.isolation.jsonl
  grep -q '"contended": false' "$iso" || continue
  jq -r --arg s "$(slug "$size")" '.scan_probes as $p | ($p.counts.summands) as $n | ($p.tsc_hz) as $hz |
    [$s] + ([$p.stages.subtract, $p.stages.key, $p.stages.filter, $p.stages.admitted, $p.stages.trial]
    | map(. / $hz * 1e9 / $n)) + [$p.counts.admitted / $n] | map(tostring) | join("\t")' "$f"
done | awk -F'\t' '{k=$1; n[k]++; for (i=2;i<=7;i++) a[k,i]+=$i}
  END {
    for (k in n) {
      s=0; for (i=2;i<=5;i++) s+=a[k,i]/n[k];
      printf "%s\t%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.4f\n", k, n[k], a[k,2]/n[k], a[k,3]/n[k], a[k,4]/n[k], a[k,5]/n[k], s, a[k,6]/n[k], a[k,7]/n[k]
    }
  }' | sort

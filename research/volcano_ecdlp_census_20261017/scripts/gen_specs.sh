#!/bin/sh
# Generates specs/*.json and runs/manifest.json (Phase B inputs). Rounds per level: ceil(2000 / workloads per round summed over D, Ua, Ub, S, N at that level).
. "$(dirname "$0")/env.sh"
cd $STUDY; mkdir -p specs runs; rm -f specs/*.json; : > $SP/t/manifest.lines
ICF=$(jq -c '.ic|{factor_base,oracle}' FROZEN.json)
for V in V3 V1 V2; do
 T=census/$V/transport
 # curves per level summed over the five measured families (N reuses D's curves)
 LC=$(for f in D floor-a floor-b crater-b D; do [ -f $T/$f.jsonl ] && cat $T/$f.jsonl; done | jq -s -c 'group_by(.level)|map({level:.[0].level,n:length})')
 for FAM in D Ua Ub S N ctl; do
  case $FAM in D) F=D;; Ua) F=floor-a;; Ub) F=floor-b;; S) F=crater-b;; N) F=D;; ctl) F=D;; esac
  [ -f $T/$F.jsonl ] || continue
  TS=20261111; [ $FAM = N ] && TS=20261115
  jq -s -c -L scripts --arg V $V --arg FAM $FAM --argjson ts $TS --argjson ic "$ICF" --argjson maxpass 4500 --argjson lvlcount "$LC" -f scripts/gen_specs.jq $T/$F.jsonl > $SP/t/gen.lines
  while IFS= read -r line; do
    n=$(echo "$line" | jq -r .name)
    echo "$line" | jq '.spec' > specs/$n.json
    echo "$line" | jq -c 'del(.spec)' >> $SP/t/manifest.lines
  done < $SP/t/gen.lines
 done
done
ls specs | wc -l

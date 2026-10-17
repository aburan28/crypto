#!/bin/sh
# Transported-instance identity: ecbench's i-th planted_shared target on each transported curve == ca's carried [k_i]G_src (hk_i).
. "$(dirname "$0")/env.sh"; cd $STUDY; mkdir -p checks; : > checks/identity.jsonl
for spec in specs/*-D-L*.json specs/*-Ua-L*.json specs/*-Ub-L*.json specs/*-S-L*.json; do
  n=$(basename $spec .json); V=${n%%-*}; fam=$(echo $n | cut -d- -f2)
  case $fam in D) F=D;; Ua) F=floor-a;; Ub) F=floor-b;; S) F=crater-b;; esac
  $ECB plan --spec $spec --json > $SP/t/plan_$n.json
  jq -c --arg n $n --arg fam $fam --slurpfile tr census/$V/transport/$F.jsonl -f scripts/identity.jq $SP/t/plan_$n.json >> checks/identity.jsonl
done
jq -s -c '{specs:length, curves:(map(.curves)|add), targets:(map(.targets)|add), mismatches:(map(.mismatches)|add)} ' checks/identity.jsonl

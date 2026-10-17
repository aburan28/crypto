#!/bin/sh
# Generate calibration/b<B>/{root.json,spec.json}: one V1-shaped root per bit size.
. "$(dirname "$0")/env.sh"
cd "$STUDY/calibration"
for B in 18 20 22 24 26; do
  mkdir -p b$B
  MIN=$((1 << (B-8)))
  line=$($VR --bits $B --ell 2 --depth 4 --seed 20261017 --min-n $MIN --trials 1 --per-level 1 --controls 0 --out $SP/t/cal$B 2>&1 | grep '^curve:')
  echo "$line" > b$B/search.txt
  P=$(echo "$line" | sed 's/.*p=\([0-9]*\) .*/\1/'); A=$(echo "$line" | sed 's/.* a=\([0-9]*\) .*/\1/'); BB=$(echo "$line" | sed 's/.* b=\([0-9]*\) .*/\1/')
  $CA volcano --p $P --a $A --b $BB --ell 2 --seed 20261017 --trials 1 --per-level 0 --controls 0 --radius 40 --max-vertices 200000 > b$B/root.json
  jq --arg b $B '. as $r | {schema:"ecbench.spec/v1", label:("volcano census calibration b="+$b), description:"Phase A calibration, one V1-shaped root, PLAN 4 A.2; known-answer synthetic curve, no security claim.",
   workloads:{curves:[{kind:"prime_explicit",p:.p,a:.curves[0].a,b:.curves[0].b,group_order:.order,r:.n,gx:(.g|split(",")[0]|tonumber),gy:(.g|split(",")[1]|tonumber)}],targets_per_curve:8,target_seed:20261017,target_kind:"planted_shared"}}' b$B/root.json > b$B/head.json
done

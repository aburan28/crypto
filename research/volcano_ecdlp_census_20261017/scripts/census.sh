#!/bin/sh
# usage: census.sh V1|V2|V3 ; writes census/<V>/walks/<src>/..., census/<V>/sources.json
. "$(dirname "$0")/env.sh"
V=$1; D=$STUDY/census/$V; SEED=20261017; TSEED=20261111; NT=8
case $V in V1) ELL=2;; V2) ELL=3;; V3) ELL=2;; esac
line=$(cat $D/root.txt)
P=$(echo "$line" | sed 's/.*p=\([0-9]*\) .*/\1/'); A=$(echo "$line" | sed 's/.* a=\([0-9]*\) .*/\1/'); B=$(echo "$line" | sed 's/.* b=\([0-9]*\) .*/\1/')
mkdir -p $D/walks/hit $D/walks/D
CAF="--ell $ELL --seed $SEED --trials 1 --per-level 0 --controls 0 --radius 40 --max-vertices 20000"
$CA volcano --p $P --a $A --b $B $CAF --census-out $D/walks/hit/vertices.jsonl > $D/walks/hit/R.json 2>$D/walks/hit/stderr.txt
CA_A=$(jq -s -r 'map(select(.level==0))|sort_by(.j)|.[0].a' $D/walks/hit/vertices.jsonl)
CA_B=$(jq -s -r 'map(select(.level==0))|sort_by(.j)|.[0].b' $D/walks/hit/vertices.jsonl)
# the scalars ecbench will draw (planted_shared depends on r and the target seed only)
jq --argjson ts $TSEED --argjson nt $NT '{schema:"ecbench.spec/v1",label:"k probe",description:"reads the planted_shared scalars",workloads:{curves:[{kind:"prime_explicit",p:.p,a:.a,b:.b,group_order:.order,r:.n,gx:(.g|split(",")[0]|tonumber),gy:(.g|split(",")[1]|tonumber)}],targets_per_curve:$nt,target_seed:$ts,target_kind:"planted_shared"},arms:[{name:"bsgs",role:"reference",method:{id:"bsgs.negation"}}],measurement:{rounds:1,warmup:0,order:"alternate",seed:1,isolation_required:"L0",timeout_seconds:60}}' $D/walks/hit/R.json > $SP/t/probe_$V.json
KS=$($ECB plan --spec $SP/t/probe_$V.json --json | jq -r '[.workloads[].planted]|map(tostring)|join(",")')
echo "$KS" > $D/ks.txt
$CA volcano --p $P --a $CA_A --b $CA_B $CAF --ks $KS --census-out $D/walks/D/vertices.jsonl --edges-out $D/walks/D/edges.jsonl --diagram $D/walks/D/diagram.json > $D/walks/D/R.json 2>$D/walks/D/stderr.txt
# choose the other sources from the D walk
H(){ printf '%s' "$1" | sha256sum | cut -c1-8 | while read h; do echo $((0x$h)); done; }
jq -s -c --argjson ha $(H "floor-a:$V") --argjson hb $(H "floor-b:$V") '
  . as $v | (map(.level)|max) as $depth |
  (map(select(.level==0))) as $cr |
  ($v | map(select(.level==$depth)) | sort_by(.j)) as $fl |
  (map(.level)|group_by(.)|map({l:.[0],c:length})|map(select(.c>=2))|.[0].l) as $ls |
  def anc($i): if $v[$i].level <= $ls then $i else anc($v[$i].parent) end;
  ($fl[$ha % ($fl|length)]) as $fa |
  ([range(0; $fl|length)] | map($fl[($hb + .) % ($fl|length)]) | map(select(anc(.id) != anc($fa.id)))|.[0]) as $fb |
  ($cr | map(select(.id != 0)) | sort_by([-.dist, .j]) | .[0]) as $cb |
  {depth:$depth, anc_level:$ls, crater_root:{id:0,a:$v[0].a,b:$v[0].b,j:$v[0].j}, floor_a:{id:$fa.id,a:$fa.a,b:$fa.b,j:$fa.j}, floor_b:{id:$fb.id,a:$fb.a,b:$fb.b,j:$fb.j},
   crater_b:(if $cb then {id:$cb.id,a:$cb.a,b:$cb.b,j:$cb.j} else null end)}' $D/walks/D/vertices.jsonl > $D/sources.json
for S in floor_a floor_b crater_b; do
  SA=$(jq -r ".$S.a // empty" $D/sources.json); [ -z "$SA" ] && continue; SB=$(jq -r ".$S.b" $D/sources.json)
  N=$(echo $S | tr _ -); mkdir -p $D/walks/$N
  $CA volcano --p $P --a $SA --b $SB $CAF --ks $KS --census-out $D/walks/$N/vertices.jsonl --edges-out $D/walks/$N/edges.jsonl > $D/walks/$N/R.json 2>$D/walks/$N/stderr.txt
done
echo "$V done"; cat $D/sources.json

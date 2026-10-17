#!/bin/sh
# Generates specs/<V>-<fam>-L<level>-c<chunk>.json and runs/manifest.json (Phase B inputs).
. "$(dirname "$0")/env.sh"
cd $STUDY; mkdir -p specs runs; rm -f specs/*.json; : > $SP/t/manifest.lines
ICF=$(jq -c '.ic' $STUDY/FROZEN.json)       # {"factor_base":..,"oracle":..,"id":..}
MAXPASS=4500
for V in V3 V1 V2; do
 T=census/$V/transport
 for FAM in D Ua Ub S N ctl; do
  case $FAM in D) F=D;; Ua) F=floor-a;; Ub) F=floor-b;; S) F=crater-b;; N) F=D;; ctl) F=D;; esac
  [ -f $T/$F.jsonl ] || continue
  TS=20261111; [ $FAM = N ] && TS=20261115
  # rounds per level from the number of workloads of the five families at that level
  jq -s -c --arg V $V --arg FAM $FAM --argjson ts $TS --argjson ic "$ICF" --argjson maxpass $MAXPASS \
     --slurpfile all <(cat $T/*.jsonl | jq -c 'select(.src!="hit")') '
   def rounds($lvl): ([$all[]|select(.level==$lvl)]|length) as $n | (2000 / ($n*8) | ceil);   # all families at that level (hit excluded in D/Ua/Ub/S; N shares D curves, counted once more below)
   . as $rows |
   ($rows | map(.level) | unique) as $levels |
   $levels[] as $L |
   ($rows | map(select(.level==$L))) as $cs |
   ([ $all[] | select(.level==$L and .src!="hit") ] | length) as $nsrc |
   (($nsrc + ([$all[]|select(.level==$L and .src=="D")]|length)) * 8) as $wlround |   # N counted as one more copy of D
   ((2000 / $wlround) | ceil) as $rounds |
   ($cs | map(
      if $FAM=="N" then {a,b,gx:(.gn|split(",")[0]|tonumber),gy:(.gn|split(",")[1]|tonumber)}
      else {a,b,gx:(.g|split(",")[0]|tonumber),gy:(.g|split(",")[1]|tonumber)} end
      | . as $c | $cs[0] | {kind:"prime_explicit",p,a:$c.a,b:$c.b,group_order:.order,r:.r,gx:$c.gx,gy:$c.gy}
      )) as $curves0 |
   (if $FAM=="ctl" then ($cs | map(. as $c | .controls[] | {kind:"prime_explicit",p:$c.p,a:.a,b:.b,group_order:$c.order,r:$c.r,gx:(.g|split(",")[0]|tonumber),gy:(.g|split(",")[1]|tonumber)})) else $curves0 end) as $curves |
   (($maxpass / (($rounds+1)*8)) | floor) as $per |
   [range(0; ($curves|length); $per)] as $starts |
   $starts | to_entries[] | . as $e |
   ($curves[$e.value : $e.value+$per]) as $chunk |
   ( [ {name:"rho-neg",role:"reference",method:{id:"rho.negation"}},
       {name:"rho-plain",role:"baseline",method:{id:"rho.plain"}},
       {name:"bsgs",role:"candidate",method:{id:"bsgs.textbook"}},
       {name:"bsgs-neg",role:"candidate",method:{id:"bsgs.negation"}},
       {name:"ic",role:"candidate",method:{id:"ic.pipeline",params:{factor_base:$ic.factor_base,oracle:$ic.oracle}}},
       {name:"rho-neg-aa",role:"control",method:{id:"rho.negation"}} ]
     + (if $V=="V3" and $L==0 and ($FAM|IN("D","N","ctl")) then [{name:"rho-cubic",role:"control",method:{id:"rho.cubic"}}] else [] end) ) as $arms |
   {name:("\($V)-\($FAM)-L\($L)-c\($e.key)"), V:$V, family:$FAM, level:$L, chunk:$e.key, curves:($chunk|length), targets:8, rounds:$rounds, warmup:1, arms:($arms|length),
    passes:(($chunk|length)*8*($rounds+1)),
    spec:{schema:"ecbench.spec/v1", label:("volcano census \($V) \($FAM) level \($L) chunk \($e.key)"),
          description:"Volcano ECDLP census Phase B (research/volcano_ecdlp_census_20261017/PLAN.md). Known-answer synthetic curves, planted_shared targets; no security claim.",
          workloads:{curves:$chunk,targets_per_curve:8,target_seed:$ts,target_kind:"planted_shared"},
          arms:$arms, measurement:{rounds:$rounds,warmup:1,order:"alternate",seed:20261017,isolation_required:"L2",timeout_seconds:600}}}
  ' $T/$F.jsonl > $SP/t/gen.lines
  while IFS= read -r line; do
    n=$(echo "$line" | jq -r .name)
    echo "$line" | jq '.spec' > specs/$n.json
    echo "$line" | jq -c 'del(.spec)' >> $SP/t/manifest.lines
  done < $SP/t/gen.lines
 done
done

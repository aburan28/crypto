#!/bin/sh
# usage: search.sh V1|V2|V3  ; logs every seed to census/search_log.jsonl; writes census/<V>/root.txt
. "$(dirname "$0")/env.sh"
V=$1
case $V in
 V1) ELL=2; DEPTH=4; BASE=20261101; EXTRA="--crater split";;
 V2) ELL=3; DEPTH=3; BASE=20261201; EXTRA="--crater split";;
 V3) ELL=2; DEPTH=4; BASE=20261301; EXTRA="--j0";;
esac
mkdir -p "$STUDY/census/$V"; LOG=$STUDY/census/search_log.jsonl
i=0
while [ $i -lt 400 ]; do
  S=$((BASE+i)); i=$((i+1))
  line=$($VR --bits 26 --ell $ELL --depth $DEPTH $EXTRA --seed $S --min-n 262144 --trials 1 --per-level 1 --controls 0 --out $SP/t/s_$V 2>&1 | grep '^curve:')
  [ -z "$line" ] && { echo "{\"V\":\"$V\",\"seed\":$S,\"result\":\"no-curve\"}" >> $LOG; continue; }
  P=$(echo "$line" | sed 's/.*p=\([0-9]*\) .*/\1/'); A=$(echo "$line" | sed 's/.* a=\([0-9]*\) .*/\1/'); B=$(echo "$line" | sed 's/.* b=\([0-9]*\) .*/\1/')
  $CA volcano --p $P --a $A --b $B --ell $ELL --seed 20261017 --trials 1 --per-level 0 --controls 0 --radius 40 --max-vertices 20000 > $SP/t/w_$V.json 2>/dev/null || { echo "{\"V\":\"$V\",\"seed\":$S,\"p\":$P,\"result\":\"ca-failed\"}" >> $LOG; continue; }
  R=$(jq -c '{p,a,b,order,n,cofactor,trace,depth,crater,vertices:.walk.vertices,per_level:.walk.per_level,special:.walk.special_j_skipped,verified:.walk.edges_verified}' $SP/t/w_$V.json)
  ok=$(jq -r --arg V $V --argjson ell $ELL --argjson P $P '
     (.walk.per_level|add) as $lv | (.walk.vertices==$lv) as $complete |
     (.cofactor<=64 and .n>=262144 and (.n % $ell)!=0 and $complete and .walk.vertices<20000) as $c |
     if $V=="V3" then ($c and .depth==4 and (($P%24)==1) and (((.n-1)%3)==0)) 
     else ($c and .crater=="split" and .walk.vertices<=400 and .walk.special_j_skipped==0 and (if $V=="V1" then .depth==4 and ($P%8)==1 else .depth==3 end)) end' $SP/t/w_$V.json)
  echo "{\"V\":\"$V\",\"seed\":$S,\"accepted\":$ok,\"walk\":$R}" >> $LOG
  if [ "$ok" = true ]; then echo "$line" > $STUDY/census/$V/root.txt; echo "$S" > $STUDY/census/$V/seed.txt; echo "$V accepted seed $S"; break; fi
done

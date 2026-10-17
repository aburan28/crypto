#!/bin/bash
# Phase B driver: resumable, manifest order. Never deletes or edits a session.
set -u
STUDY=/home/user/crypto/research/volcano_ecdlp_census_20261017
ECB=/home/user/crypto/target/release/ecbench
cd "$STUDY" || exit 2
PROG=runs/progress.jsonl
MAXWALL=$((16*3600)); MINFREE_KB=$((1536*1024)); MAXFAIL=3
START=$(jq -r 'select(.type=="header").start_unix' $PROG | head -1)
now() { date -u +%Y-%m-%dT%H:%M:%SZ; }
done_p() { jq -e --arg s "$1" 'select(.spec==$s and .status=="done")' $PROG >/dev/null 2>&1 && return 0; return 1; }
free_kb() { df -Pk "$STUDY" | awk 'NR==2{print $4}'; }
audit_streak() { jq -r 'select(.spec!=null)|.audit_ok' $PROG | tac | awk 'BEGIN{n=0}{if($1=="true")exit; n++}END{print n}'; }
log() { jq -nc --arg t "$(now)" --arg m "$1" '{type:"note",time:$t,msg:$m}' >> $PROG; }

for name in $(jq -r '.sessions[].name' runs/manifest.json); do
  done_p "$name" && continue
  el=$(( $(date +%s) - START ))
  if [ $el -ge $MAXWALL ]; then log "guard: 16h wall reached before $name"; echo GUARD_WALL; exit 3; fi
  if [ "$(free_kb)" -lt $MINFREE_KB ]; then log "guard: free disk < 1.5 GiB before $name"; echo GUARD_DISK; exit 3; fi
  if [ "$(audit_streak)" -ge $MAXFAIL ]; then log "guard: 3 failed audits in a row before $name"; echo GUARD_AUDIT; exit 3; fi
  out="sessions/$name"; n=1
  # an existing directory is an interrupted attempt: keep it, record it, rerun to .rN
  while [ -e "$out" ]; do
    if ! jq -e --arg o "$out" 'select(.out==$o)' $PROG >/dev/null 2>&1; then
      jq -nc --arg t "$(now)" --arg s "$name" --arg o "$out" '{type:"orphan_dir",time:$t,spec:$s,out:$o,note:"directory existed without a progress line; kept as evidence"}' >> $PROG
    fi
    n=$((n+1)); out="sessions/$name.r$n"
  done
  t0=$(date +%s)
  $ECB run --wait --cpus auto --spec specs/$name.json --out $out > $out.log 2>&1; rc_run=$?
  # ecbench creates the dir; the log was opened before it. Keep both.
  $ECB verify --dir $out --replay 12 --out $out.audit.json > $out.verify.log 2>&1; rc_ver=$?
  audit_ok=$(jq -r '.ok' $out.audit.json 2>/dev/null || echo false)
  $ECB table --dir $out > runs/tables/$(basename $out).txt 2>&1
  base=$(basename $out)
  # compact per-run extract
  { echo 'spec,volcano,family,level,chunk,curve_slug,workload_id,target_index,arm,role,round,warmup,ops_gae,S,online_wall_ns,solve_wall_ns,status,verified,isolation'
    jq -r --arg spec "$base" --arg V "$(jq -r --arg n $name '.sessions[]|select(.name==$n)|.V' runs/manifest.json)" --arg F "$(jq -r --arg n $name '.sessions[]|select(.name==$n)|.family' runs/manifest.json)" --arg L "$(jq -r --arg n $name '.sessions[]|select(.name==$n)|.level' runs/manifest.json)" --arg C "$(jq -r --arg n $name '.sessions[]|select(.name==$n)|.chunk' runs/manifest.json)" \
      '[$spec,$V,$F,$L,$C,.workload.curve.slug,.workload.workload_id,.workload.target_index,.arm,.role,.round,.warmup,.cost.total_gae,.cost.s,.online.wall_ns,.time.solve_wall_ns,.outcome.status,(.outcome.status=="verified"),.isolation.level]|@csv' $out/records.jsonl
  } | xz -9e > runs/extract/$base.csv.xz
  runs=$(wc -l < $out/records.jsonl)
  fails=$(jq -c 'select(.warmup==false and .outcome.status!="verified")' $out/records.jsonl | wc -l)
  sstat=$(jq -r '.status' $out/session.json)
  raw=$(du -sb $out | cut -f1)
  tar -C sessions -cf - $base | xz -9e -T1 > sessions/$base.tar.xz
  xzb=$(stat -c %s sessions/$base.tar.xz); sha=$(sha256sum sessions/$base.tar.xz | cut -d' ' -f1)
  wall=$(( $(date +%s) - t0 ))
  st=done; [ "$sstat" != "complete" ] && st=$sstat; [ "$audit_ok" != "true" ] && st=audit_failed
  [ "$sstat" != "complete" ] && audit_ok=false
  jq -nc --arg spec "$name" --arg out "$out" --arg st "$st" --arg ss "$sstat" --argjson ok "$audit_ok" --argjson runs $runs --argjson fails $fails \
    --argjson raw $raw --argjson xz $xzb --arg sha "$sha" --argjson wall $wall --arg t "$(now)" --argjson rr $rc_run --argjson rv $rc_ver \
    '{spec:$spec,out:$out,status:$st,session_status:$ss,audit_ok:$ok,runs:$runs,failures:$fails,bytes_raw:$raw,bytes_xz:$xz,sha256_tar_xz:$sha,wall_s:$wall,time:$t,rc_run:$rr,rc_verify:$rv}' >> $PROG
done
log "manifest finished"
echo ALLDONE

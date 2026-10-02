#!/usr/bin/env bash
# Recompute the sigma inline A/B result from a fetched raw result directory.
set -uo pipefail

R=${1:?usage: audit.sh RESULT_DIR [--allow-incomplete]}
allow_incomplete=${2:-}
fail=0

hash_file() {
  if command -v sha256sum >/dev/null 2>&1; then
    sha256sum "$1" | awk '{print $1}'
  else
    shasum -a 256 "$1" | awk '{print $1}'
  fi
}

need_file() {
  if [ ! -s "$R/$1" ]; then
    echo "missing or empty: $1" >&2
    fail=1
  fi
}

for path in config.txt host-before.txt source-files.sha256 source-manifest.sha256 \
  binary-sha256.txt code-sha256.txt resources-control.txt resources-inline3.txt \
  sass-control.txt.gz sass-inline3.txt.gz gate-control-arithmetic.log \
  gate-inline3-arithmetic.log gate-control-storage.log gate-inline3-storage.log \
  gate-control-shared-sigma.log gate-inline3-shared-sigma.log checkpoint-sha256.txt \
  verify-control.log verify-inline3.log dp-control.bin dp-inline3.bin dp-identity.txt; do
  need_file "$path"
done

for name in control inline3; do
  grep -qx 'PASS: 6240 GPU paired Frobenius vectors, both inputs against independent routing' \
    "$R/gate-$name-arithmetic.log" 2>/dev/null || fail=1
  grep -qx 'PASS: 128 GPU storage cases, 297344 records, independent physical images and logical reads with canaries' \
    "$R/gate-$name-storage.log" 2>/dev/null || fail=1
  grep -qx 'PASS: 114 complete block mask snapshots, 51072 words, output guards and inactive blocks' \
    "$R/gate-$name-shared-sigma.log" 2>/dev/null || fail=1
  grep -Eq '\(300 verified against the reference, 0 dropped\)$' \
    "$R/verify-$name.log" 2>/dev/null || fail=1
  if grep -Eq 'MISMATCH|OVERFLOW|stopping:|collision found|verified \[k\]P|solved' \
      "$R/verify-$name.log" 2>/dev/null; then fail=1; fi
done
grep -Eq '^inline3 +[1-9][0-9]* records .* IDENTICAL to control$' \
  "$R/dp-identity.txt" 2>/dev/null || fail=1

if [ ! -s "$R/samples.tsv" ]; then
  fail=1
  if [ "$allow_incomplete" = '--allow-incomplete' ]; then
    printf '{\n  "valid": false,\n  "decision": "no-go",\n  "reason": "pre-timing gate failed; timing suppressed"\n}\n' > "$R/result.json"
    exit 0
  fi
fi

if [ "$fail" -eq 0 ]; then
  warmups=$(awk -F '\t' 'NR>1 && $1=="warmup" {n++} END {print n+0}' "$R/samples.tsv")
  ranked=$(awk -F '\t' 'NR>1 && $1=="ranked" {n++} END {print n+0}' "$R/samples.tsv")
  [ "$warmups" -eq 2 ] || fail=1
  [ "$ranked" -eq 10 ] || fail=1
  awk -F '\t' 'NR>1 && ($6 != 201863462912 || $7 != 0 || $5+0 <= 0) {bad=1} END {exit bad}' \
    "$R/samples.tsv" || fail=1
  # Require exactly two arms in every pair and the frozen alternating order.
  awk -F '\t' '
    NR>1 && $1=="ranked" {
      seen[$2,$4]++; order[$2,$3]=$4; count[$2]++
    }
    END {
      for (p=1;p<=5;p++) {
        if (count[p]!=2 || seen[p,"control"]!=1 || seen[p,"inline3"]!=1) bad=1
        want1=(p%2 ? "control" : "inline3"); want2=(p%2 ? "inline3" : "control")
        if (order[p,1]!=want1 || order[p,2]!=want2) bad=1
      }
      exit bad
    }' "$R/samples.tsv" || fail=1
fi

if [ "$fail" -eq 0 ]; then
  awk -F '\t' '
    NR>1 && $1=="ranked" {
      rate[$2,$4]=$5+0
    }
    END {
      for (p=1;p<=5;p++) {
        c[p]=rate[p,"control"]; x[p]=rate[p,"inline3"]; r[p]=x[p]/c[p]
        lc[p]=log(r[p]); sumc+=c[p]; sumx+=x[p]; suml+=lc[p]
      }
      for (i=1;i<=5;i++) for (j=i+1;j<=5;j++) {
        if (c[j]<c[i]) {t=c[i];c[i]=c[j];c[j]=t}
        if (x[j]<x[i]) {t=x[i];x[i]=x[j];x[j]=t}
        if (r[j]<r[i]) {t=r[i];r[i]=r[j];r[j]=t}
      }
      meanl=suml/5; ss=0
      for (p=1;p<=5;p++) ss+=(lc[p]-meanl)^2
      sd=sqrt(ss/4); half=2.7764451051977987*sd/sqrt(5)
      printf "controlMedianM=%.9f\n",c[3]
      printf "candidateMedianM=%.9f\n",x[3]
      printf "medianRatio=%.12f\n",r[3]
      printf "minimumRatio=%.12f\n",r[1]
      printf "maximumRatio=%.12f\n",r[5]
      printf "geometricMeanRatio=%.12f\n",exp(meanl)
      printf "ciLower=%.12f\n",exp(meanl-half)
      printf "ciUpper=%.12f\n",exp(meanl+half)
      printf "candidateRatioTo26B=%.12f\n",x[3]/26000.0
      printf "allPairsWin=%s\n",(r[1]>1 ? "true" : "false")
    }' "$R/samples.tsv" > "$R/statistics.env"
  # The file contains only fixed names and numeric/boolean values emitted above.
  # shellcheck disable=SC1090,SC1091
  . "$R/statistics.env"
  if [ "$allPairsWin" = true ] && awk -v lo="$ciLower" 'BEGIN {exit !(lo>1)}'; then
    decision=promote
  else
    decision=no-go
  fi
else
  controlMedianM=0 candidateMedianM=0 medianRatio=0 minimumRatio=0 maximumRatio=0
  geometricMeanRatio=0 ciLower=0 ciUpper=0 candidateRatioTo26B=0 allPairsWin=false
  decision=no-go
fi

source_rev=$(sed -n 's/^source=//p' "$R/config.txt" 2>/dev/null | head -1)
gpu_line=$(head -1 "$R/host-before.txt" 2>/dev/null)
source_manifest=$(hash_file "$R/source-files.sha256" 2>/dev/null || true)
binary_manifest=$(hash_file "$R/binary-sha256.txt" 2>/dev/null || true)
code_manifest=$(hash_file "$R/code-sha256.txt" 2>/dev/null || true)
dp_identity=$(hash_file "$R/dp-identity.txt" 2>/dev/null || true)
samples_digest=$(hash_file "$R/samples.tsv" 2>/dev/null || true)

if ! command -v jq >/dev/null 2>&1; then
  echo 'jq is required to emit result.json' >&2
  exit 1
fi

samples_json=$(jq -Rn '[inputs | select(length>0) | split("\t") |
  select(.[0] != "phase") |
  {phase:.[0], pair:(.[1]|tonumber), order:(.[2]|tonumber), variant:.[3],
   rateMps:(.[4]|tonumber), iterations:(.[5]|tonumber), dropped:(.[6]|tonumber),
   logSha256:.[7], gpuState:.[8]}]' < "$R/samples.tsv")

jq -n \
  --argjson valid "$([ "$fail" -eq 0 ] && echo true || echo false)" \
  --arg decision "$decision" \
  --arg sourceRev "$source_rev" \
  --arg gpu "$gpu_line" \
  --arg sourceManifestSha256 "$source_manifest" \
  --arg binaryManifestSha256 "$binary_manifest" \
  --arg codeManifestSha256 "$code_manifest" \
  --arg dpIdentitySha256 "$dp_identity" \
  --arg samplesSha256 "$samples_digest" \
  --argjson samples "$samples_json" \
  --argjson controlMedianM "$controlMedianM" \
  --argjson candidateMedianM "$candidateMedianM" \
  --argjson medianRatio "$medianRatio" \
  --argjson minimumRatio "$minimumRatio" \
  --argjson maximumRatio "$maximumRatio" \
  --argjson geometricMeanRatio "$geometricMeanRatio" \
  --argjson ciLower "$ciLower" \
  --argjson ciUpper "$ciUpper" \
  --argjson candidateRatioTo26B "$candidateRatioTo26B" \
  --argjson allPairsWin "$allPairsWin" \
  '{schema:"ecc2k130-sigma-inline-refresh-v1", valid:$valid,
    class:"engineering stage diagnostic", decision:$decision,
    claimBoundary:"Complete scalar sigma-walk throughput on one Modal RTX PRO 6000; no search, solver or full-ECDLP speedup claim",
    modes:{control:{packedInlinePoly:0},candidate:{name:"inline3",packedInlinePoly:3}},
    fixed:{batch:16,blockThreads:256,minBlocks:2,stateTile:256,compactState:1,
      sharedSigma:1,nativeCarrylessMultiply:1,workers:385024,scalarWalks:6160384,
      runId:7,steps:1024,launches:32,updatesPerRankedRow:201863462912,
      compiler:"CUDA 13.3.73 sm_120",warmupsExcluded:2,alternatingPairs:5},
    correctness:{arithmetic:true,storage:true,sharedSigma:true,
      checkpointCrossBinary:true,referenceReplaysPerArm:300,
      replayMismatchCount:0,droppedCount:0,formatAwareSortedDpIdentity:true,
      deterministicSpreadReplayAvailable:false},
    timing:{unit:"million complete scalar updates per second",
      controlMedianMps:$controlMedianM,candidateMedianMps:$candidateMedianM,
      medianPairedRatio:$medianRatio,minimumPairedRatio:$minimumRatio,
      maximumPairedRatio:$maximumRatio,geometricMeanPairedRatio:$geometricMeanRatio,
      pairedLogRatio95CI:[$ciLower,$ciUpper],allPairsFavorCandidate:$allPairsWin,
      candidateRatioTo26BObjective:$candidateRatioTo26B,
      objective26BAchieved:($candidateMedianM >= 26000)},
    admission:{requiresAllPairsAboveOne:true,requiresPaired95CILowerAboveOne:true},
    provenance:{sourceRev:$sourceRev,gpu:$gpu,sourceManifestSha256:$sourceManifestSha256,
      binaryManifestSha256:$binaryManifestSha256,codeManifestSha256:$codeManifestSha256,
      dpIdentitySha256:$dpIdentitySha256,samplesSha256:$samplesSha256,
      rawFiles:{source:"source-files.sha256",binaries:"binary-sha256.txt",
        code:"code-sha256.txt",resources:["resources-control.txt","resources-inline3.txt"],
        sass:["sass-control.txt.gz","sass-inline3.txt.gz"],samples:"samples.tsv"}},
    samples:$samples}' > "$R/result.json"

if [ "$fail" -ne 0 ] && [ "$allow_incomplete" != '--allow-incomplete' ]; then
  exit 1
fi

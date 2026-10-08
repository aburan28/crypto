#!/bin/bash
# n=61 L=65536 compact-orbit vs frozen KS v2 batched rho: restart of the
# aborted growing_n_n61_L65536.sh panel. Resumable: every stage writes a
# .done marker and is skipped when rerun, so an abort keeps the corpus and
# finished blocks.
#
# Deviations from growing_n_n61_L65536.sh (recorded in growing_n_n61_L65536_20260930/RESULT.md):
# - corpora come from KIC_RHO_GENERATE_ONLY=1 (same blake3 fixture derivation
#   as the walking run) instead of a full rho walk per corpus;
# - K candidates are skipped when the estimated IC peak RSS (94 B per regular
#   state, K^2*n states, from the L=16384 K=1400 run) exceeds available memory;
# - tools/isolated_bench.py does not run on macOS (no sched_getaffinity), so
#   walls come from /usr/bin/time -l with load and memory recorded per run and
#   are not AGENTS.md section-10 evidence.
# No L=32 panels.
set -u
LAB=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$LAB/../../.." && pwd)
OUT=${OUT_OVERRIDE:-$LAB/growing_n_n61_L65536_20260930}
mkdir -p "$OUT"
cd "$OUT"
KS=${KS_OVERRIDE:-$ROOT/target/release/examples/koblitz_rho_batch_ks_v2_n61}
IC=${IC_OVERRIDE:-$ROOT/target/release/examples/koblitz_orbit_dlp_fast}
L=65536
N=61
SEED=531310
BLOCKS=${BLOCKS:-3}
candidates=${KS_CANDIDATES:-"1400 1800 2200"}

log() { echo "$* $(date -u +%FT%TZ)" | tee -a panel.log; }
wall() { grep ' real' "$1" | awk '{print $1}'; }
avail_bytes() {
  vm_stat | awk '/page size/{ps=$8} /Pages free/{f=$3} /Pages inactive/{i=$3}
    /Pages purgeable/{p=$3} /Pages speculative/{s=$3}
    END{gsub(/\./,"",f);gsub(/\./,"",i);gsub(/\./,"",p);gsub(/\./,"",s);
        printf "%d\n",(f+i+p+s)*ps}'
}
load_now() { sysctl -n vm.loadavg | tr -d '{}' | awk '{print $1, $2, $3}'; }
fits() {
  local K=$1 need avail
  need=$((94 * N * K * K))
  avail=$(avail_bytes)
  echo "K=$K est_rss=$need avail=$avail" >> memory_checks.log
  [ "$need" -lt "$avail" ]
}

if [ ! -s host.json ]; then
  python3 - "$ROOT" > host.json <<'PY'
import json, platform, subprocess, sys
def sh(*a):
    return subprocess.run(a, capture_output=True, text=True, cwd=sys.argv[1]).stdout.strip()
print(json.dumps({
    "git_head": sh("git", "rev-parse", "HEAD"),
    "rustc": sh("rustc", "--version"),
    "cpu": sh("sysctl", "-n", "machdep.cpu.brand_string"),
    "logical_cpus": int(sh("sysctl", "-n", "hw.ncpu")),
    "memory_bytes": int(sh("sysctl", "-n", "hw.memsize")),
    "os": platform.platform(), "machine": platform.machine(),
    "isolated_bench": "unsupported on macOS: os.sched_getaffinity missing",
}, indent=1))
PY
fi
shasum -a 256 "$KS" "$IC" > SHA256SUMS
log "start KS=$KS IC=$IC L=$L N=$N BLOCKS=$BLOCKS"

corpus() {
  local name=$1
  if [ ! -f "corpus_$name.done" ]; then
    KIC_RHO_GENERATE_ONLY=1 KIC_RHO_BATCH_CORPUS=$name \
      "$KS" $N 0 signed_frobenius $L $SEED > "corpus_$name.jsonl"
    python3 -c "import json; [print(d['published_fixture_scalar']) for d in map(json.loads, open('corpus_$name.jsonl')) if d.get('kind')=='rho_ks_public_fixture']" \
      > "scalars_$name.txt"
    [ "$(wc -l < "scalars_$name.txt")" -eq $L ] || { log "corpus $name short"; exit 1; }
    shasum -a 256 "scalars_$name.txt" >> SHA256SUMS.corpora
    rm -f "corpus_$name.jsonl"
    touch "corpus_$name.done"
    log "corpus $name ready"
  fi
}

tune=n61-ks-growing-tune-$L-v1
eval_name=n61-ks-growing-$L-v1
corpus "$tune"
corpus "$eval_name"

if [ ! -s chosen_K_n61.txt ]; then
  for K in $candidates; do
    [ -f "tune_n61_K$K.done" ] && continue
    if ! fits "$K"; then
      log "tune n=61 K=$K skipped: estimated RSS exceeds available memory"
      echo skipped > "tune_n61_K$K.done"
      continue
    fi
    log "tune K=$K start load=$(load_now)"
    /usr/bin/time -l "$IC" "construct:$N:0:$K" "scalars_$tune.txt" 7 /dev/null \
      > "tune_n61_K$K.summary.json" 2> "tune_n61_K$K.time"
    rc=$?
    log "tune n=61 K=$K rc=$rc wall=$(wall "tune_n61_K$K.time") load=$(load_now)"
    [ $rc = 0 ] && echo ok > "tune_n61_K$K.done"
  done
  best=$(for K in $candidates; do
    [ "$(cat "tune_n61_K$K.done" 2>/dev/null)" = ok ] && echo "$K $(wall "tune_n61_K$K.time")"
  done | sort -k2 -g | head -1 | awk '{print $1}')
  [ -n "$best" ] || { log "no K completed"; exit 1; }
  echo "$best" > chosen_K_n61.txt
  log "chosen n=61 K=$best"
fi
K=$(cat chosen_K_n61.txt)

run_ks() {
  local b=$1
  [ -f "ks_n61_b$b.done" ] && return
  load_now > "ks_n61_b$b.load"
  log "ks b=$b start load=$(cat "ks_n61_b$b.load")"
  KIC_RHO_DP_BITS=4 KIC_RHO_BATCH_CORPUS=$eval_name \
    /usr/bin/time -l "$KS" $N 0 signed_frobenius $L $SEED \
    > "ks_n61_b$b.jsonl" 2> "ks_n61_b$b.time"
  local rc=$?
  log "ks n=61 b=$b rc=$rc wall=$(wall "ks_n61_b$b.time")"
  [ $rc = 0 ] && touch "ks_n61_b$b.done"
}
run_ic() {
  local b=$1
  [ -f "ic_n61_b$b.done" ] && return
  load_now > "ic_n61_b$b.load"
  avail_bytes > "ic_n61_b$b.avail"
  log "ic b=$b K=$K start load=$(cat "ic_n61_b$b.load")"
  KIC_DUMP_BASE=$OUT/base_n61_K$K.jsonl \
    /usr/bin/time -l "$IC" "construct:$N:0:$K" "scalars_$eval_name.txt" 7 "ic_n61_b$b.jsonl" \
    > "ic_n61_b$b.summary.json" 2> "ic_n61_b$b.time"
  local rc=$?
  log "ic n=61 K=$K b=$b rc=$rc wall=$(wall "ic_n61_b$b.time")"
  [ $rc = 0 ] && touch "ic_n61_b$b.done"
}

for b in $(seq 0 $((BLOCKS - 1))); do
  if [ $((b % 2)) = 0 ]; then run_ic $b; run_ks $b; else run_ks $b; run_ic $b; fi
done
log "BLOCKS_DONE"

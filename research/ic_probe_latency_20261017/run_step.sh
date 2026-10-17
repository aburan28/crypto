#!/bin/sh
# One step of the IC probe-latency study. Frozen by PROTOCOL.md; the
# autoresearcher trial plan calls this once per trial with a fresh OUT dir.
#
#   run_step.sh STEP OUT [ARG]
#
# Steps: build, stage0, sweep, ic-4k, ic-2m, rho, claims RUNS...
# Every step refuses a non-empty OUT, checks the frozen source and input
# hashes, and records the binaries' SHA-256.
set -eu

study=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
root=$(CDPATH= cd -- "$study/../.." && pwd)
step=${1:?step}
out=${2:?output directory}
n83=$root/research/n83_dense_root_20261013
producer=$root/target/release/examples/koblitz_orbit_dlp_fast
tool=$root/target/release/examples/ic_probe_latency
rho=$root/target/release/examples/n83_rho_optimization_eval

mkdir -p "$out"
[ -z "$(ls -A "$out")" ] || { echo "output directory not empty: $out" >&2; exit 2; }

sha() { sha256sum "$1" | awk '{print $1}'; }

# Frozen sources (PROTOCOL.md "Frozen identities") and inputs.
check_sources() {
  while read -r want path; do
    got=$(sha "$root/$path")
    [ "$got" = "$want" ] || { echo "source hash mismatch: $path" >&2; exit 2; }
  done < "$study/sources.sha256"
  [ "$(sha "$n83/targets.jsonl")" = e9ed7ddb0aeafdbfc7354d67279a8712572eea9d2aee3896f99633bdc2e937cc ] \
    || { echo "target panel hash mismatch" >&2; exit 2; }
  [ "$(sha "$n83/answers.jsonl")" = 6a7022f739df0b2a9e7d7bbbf542820244497147f98d595dff0dbafba0bb9181 ] \
    || { echo "answer sidecar hash mismatch" >&2; exit 2; }
}

stamp() {
  date -u +%Y-%m-%dT%H:%M:%SZ > "$out/started_utc.txt"
  {
    printf 'step\t%s\n' "$step"
    printf 'commit\t%s\n' "$(git -C "$root" rev-parse HEAD)"
    for b in "$producer" "$tool" "$rho"; do
      [ -x "$b" ] && printf 'sha256\t%s\t%s\n' "$(basename "$b")" "$(sha "$b")"
    done
    printf 'thp\t%s\n' "$(cat /sys/kernel/mm/transparent_hugepage/enabled 2>/dev/null || echo none)"
    printf 'cpu\t%s\n' "$(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ //')"
    printf 'nproc\t%s\n' "$(nproc)"
    printf 'mem_total_kib\t%s\n' "$(awk '/MemTotal/ {print $2}' /proc/meminfo)"
  } > "$out/launch.tsv"
}

finish() {
  printf '%s\n' "$1" > "$out/exit_status.txt"
  date -u +%Y-%m-%dT%H:%M:%SZ > "$out/finished_utc.txt"
  echo "$step status=$1"
  exit "$1"
}

# chase BYTES PAGES FILE
chase() {
  "$tool" chase --bytes "$1" --pages "$2" --accesses 16Mi --reps 7 --seed 20261017 --out "$3"
}

# producer K TARGETS HUGE CAP
ic() {
  set +e
  timeout "$4" env KIC_PACKED_ROOT_TABLE=1 KIC_PACKED_ROOT_TABLE_DENSE=1 KIC_ROOT_BLOOM=0 \
    KIC_RANK_THREADS=4 KIC_TARGET_THREADS=1 KIC_ROOT_TABLE_HUGEPAGE="$3" \
    "$producer" "construct:83:1:$1" "$2" 29 "$out/targets_result.jsonl" \
    > "$out/stdout.jsonl" 2> "$out/stderr.log"
  status=$?
  set -e
  return $status
}

check_sources
case "$step" in
  build)
    stamp
    ( cd "$root" && cargo build --release --locked \
        --example koblitz_orbit_dlp_fast --example ic_probe_latency \
        --example n83_rho_optimization_eval ) > "$out/build.log" 2>&1 || finish $?
    rm -f "$out/launch.tsv"; stamp
    finish 0 ;;
  stage0)
    stamp
    a=$n83/runs/panel
    "$tool" geomfit \
      --label k1271_regular --targets "$a/k1271_regular/targets_result.jsonl" \
      --label k1271_dense --targets "$a/k1271_dense/targets_result.jsonl" \
      --label k1500_dense --targets "$a/k1500_dense/targets_result.jsonl" \
      --label k1798_dense --targets "$a/k1798_dense/targets_result.jsonl" \
      --out "$out/geomfit.json" 2> "$out/stderr.log" || finish $?
    finish 0 ;;
  sweep)
    stamp
    for bytes in 1Mi 16Mi 256Mi 1Gi 8Gi; do
      for pages in 4k 2m; do
        chase "$bytes" "$pages" "$out/chase-$bytes-$pages.json" 2>> "$out/stderr.log" || finish $?
      done
    done
    finish 0 ;;
  ic-4k|ic-2m)
    stamp
    huge=0; pages=4k
    [ "$step" = ic-2m ] && { huge=1; pages=2m; }
    # The null object runs first, on the same page size, in the same trial.
    chase 8Gi "$pages" "$out/chase.json" 2> "$out/chase.stderr.log" || finish $?
    ic 1798 "$n83/targets.jsonl" "$huge" 10800 || finish $?
    finish 0 ;;
  rho)
    stamp
    set +e
    timeout 600 "$rho" "$n83/primary_target.jsonl" 2026101301 "$n83/rho_config.json" \
      "$out/receipt.json" > "$out/stdout.log" 2> "$out/stderr.log"
    status=$?
    set -e
    finish $status ;;
  claims)
    # claims OUT STAGE0 A1 B1 A2 B2 RHO1 RHO2
    shift 2
    [ $# -eq 7 ] || { echo "claims needs 7 run directories" >&2; exit 2; }
    stamp
    "$tool" claims --stage0 "$1/geomfit.json" \
      --archived "$n83/runs/panel/k1798_dense/targets_result.jsonl" \
      --answers "$n83/answers.jsonl" \
      --a1 "$2" --b1 "$3" --a2 "$4" --b2 "$5" --rho "$6,$7" \
      --out "$out/claims.json" 2> "$out/stderr.log" || finish $?
    finish 0 ;;
  *)
    echo "unknown step: $step" >&2; exit 2 ;;
esac

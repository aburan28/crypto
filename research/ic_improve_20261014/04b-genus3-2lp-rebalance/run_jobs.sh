#!/usr/bin/env bash
# Thin orchestration for PROTOCOL.md §9: run each job line "<name> <tool args...>" under
# timeout 3600 s and a 3 GiB address-space cap, four at a time. A job whose output file
# already exists is skipped (resume). Failures go to <outdir>/failures.jsonl with the reason.
set -euo pipefail
BIN=${BIN:?set BIN to the compiled ic_improve_04b_genus3_2lp binary}
JOBS=$1; OUTDIR=$2; PAR=${PAR:-4}
mkdir -p "$OUTDIR"
run_one() {
  name=$1; shift
  out="$OUTDIR/$name.jsonl"; [ -s "$out" ] && return 0
  start=$(date -u +%Y-%m-%dT%H:%M:%SZ); t0=$(date +%s.%N)
  set +e
  prlimit --as=3221225472 timeout 3600 "$BIN" "$@" > "$out.tmp" 2> "$OUTDIR/$name.err"; rc=$?
  set -e
  wall=$(echo "$(date +%s.%N) - $t0" | bc)
  if [ $rc -eq 0 ] && [ -s "$out.tmp" ]; then mv "$out.tmp" "$out"; [ -s "$OUTDIR/$name.err" ] || rm -f "$OUTDIR/$name.err"
  else
    reason="exit $rc"; [ $rc -eq 124 ] && reason="timeout 3600 s"
    printf '{"job":"%s","args":"%s","exit":%d,"reason":"%s","start":"%s","wall_s":%s}\n' "$name" "$*" "$rc" "$reason" "$start" "$wall" >> "$OUTDIR/failures.jsonl"
    rm -f "$out.tmp"
  fi
}
export -f run_one; export BIN OUTDIR
xargs -P "$PAR" -L 1 bash -c 'run_one "$@"' _ < "$JOBS"

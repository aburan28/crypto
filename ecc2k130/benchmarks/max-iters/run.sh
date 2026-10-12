#!/usr/bin/env bash
# PROTOCOL.md end to end: build the scale-model client, walk the four arms on
# each run id, analyze, and write results.json, results.txt, manifest.json
# and raw/*.log here.
#
#   ./run.sh                         all three run ids, four runs at a time
#   RUN_IDS="30001" ./run.sh walk    the development run id only, no analysis
#   ./run.sh analyze                 analyze an existing WORK directory
#
# The corpora and checkpoints (about 20 MB and 1.3 MB per run) stay in WORK.
# A single-thread run is deterministic in its run id, so run.sh regenerates
# them byte for byte; manifest.json holds their SHA-256 to check that against.
set -euo pipefail
HERE=$(cd "$(dirname "$0")" && pwd)
ROOT=$(cd "$HERE/../.." && pwd)
WORK=${WORK:-${TMPDIR:-/tmp}/ecc2k130-max-iters}
JOBS=${JOBS:-4}
RUN_IDS=${RUN_IDS:-"30001 30002 30003"}
ARMS="ref c30 c31 c32"
BIN="$WORK/ecc2k130-cpu-guard64"
mkdir -p "$WORK"

cap() {
    case "$1" in ref) echo 0 ;; c30) echo 20800 ;; c31) echo 41600 ;; c32) echo 83200 ;; esac
}

build() {
    # The production host build (`make cpu`), with one macro added.
    local cmd
    cmd=$(cd "$ROOT" && make -s -n -B ecc2k130-cpu | tr '\n' ' ' | sed 's/\\ / /g')
    if [ "$(grep -o ' -o ' <<<"$cmd" | wc -l)" != 1 ] || [[ "$cmd" != *" -o ecc2k130-cpu"* ]]; then
        echo "expected one compile of ecc2k130-cpu, got: $cmd" >&2
        exit 1
    fi
    cmd=${cmd/ -o ecc2k130-cpu/ -DECC_GUARD_PERIOD=64 -o $BIN}
    echo "$cmd" > "$WORK/build.cmd"
    (cd "$ROOT" && eval "$cmd")
}

walkOne() {
    local runId=$1 arm=$2 base="$WORK/$1-$2" t0 rc
    rm -f "$base.bin" "$base.ck" "$base.ck.tmp"
    t0=$(date +%s.%N)
    set +e
    OMP_NUM_THREADS=1 "$BIN" --curve 131 --threads 1 --steps 1024 --launches 108 \
        --dp-weight 44 --dp-cap 65536 --verify 4 --run-id "$runId" \
        --max-iters "$(cap "$arm")" --dp-file "$base.bin" --checkpoint "$base.ck" \
        > "$base.log" 2>&1
    rc=$?
    set -e
    echo "$rc $(awk -v a="$t0" -v b="$(date +%s.%N)" 'BEGIN { printf "%.1f", b - a }')" > "$base.status"
    echo "$runId $arm exit $rc"
}

walk() {
    [ -x "$BIN" ] || build
    export -f walkOne cap
    export WORK BIN
    for runId in $RUN_IDS; do for arm in $ARMS; do echo "$runId $arm"; done; done |
        xargs -P "$JOBS" -L 1 bash -c 'walkOne "$0" "$1"'
}

analyze() {
    python3 "$HERE/analyze.py" --work "$WORK" --out "$HERE"
    mkdir -p "$HERE/raw"
    cp "$WORK"/*.log "$WORK"/*.status "$WORK/build.cmd" "$HERE/raw/"
    python3 - "$WORK" "$HERE" "$ROOT" "$BIN" <<'EOF'
import hashlib, json, os, subprocess, sys
work, here, root, binary = sys.argv[1:]
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()
git = lambda *a: subprocess.check_output(["git", *a], cwd=root, text=True).strip()
files = sorted(f for f in os.listdir(work) if f.endswith((".bin", ".ck")))
manifest = {
    "protocol": "PROTOCOL.md",
    "sourceCommit": git("rev-parse", "HEAD"),
    "sourceDirty": bool(git("status", "--porcelain", "--", "include", "src", "generated", "Makefile")),
    "buildCommand": open(os.path.join(work, "build.cmd")).read().strip(),
    "binarySha256": sha(binary),
    "raw": {f: {"bytes": os.path.getsize(os.path.join(work, f)), "sha256": sha(os.path.join(work, f))}
            for f in files},
    "regenerate": "WORK=<dir> ./run.sh (single-thread runs are deterministic in the run id)",
}
with open(os.path.join(here, "manifest.json"), "w") as fh:
    json.dump(manifest, fh, indent=1, sort_keys=True)
    fh.write("\n")
EOF
}

case "${1:-all}" in
    build) build ;;
    walk) walk ;;
    analyze) analyze ;;
    all) walk; analyze ;;
    *) echo "usage: $0 [all|build|walk|analyze]" >&2; exit 2 ;;
esac

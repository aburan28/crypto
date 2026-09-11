#!/usr/bin/env bash
#
# Rehearse the distributed campaign's end game on a curve that finishes in
# seconds: several workers with disjoint run ids each write their own corpus
# file, none of them holds both halves of the collision, and merge.py has to
# find the pair across files and recover k through the host client.
#
# Same idea as ../rehearse.sh, but the merge is the one the EC2 campaign uses
# (bucketed external merge) rather than the client's in-process --load, and
# the walks are shaped like the real ones: each lives several distinguished-
# point intervals (steps x launches >> 1/p) so a walk that merges into another
# trail actually reaches the shared point and reports it.  On GF(2^41) with
# cutoff 18 about a fifth of points are distinguished, so eight-step launches
# repeated a few times give walks of 3-6 intervals; the launch count is the
# budget knob.  Too much total work and one worker solves alone, too little
# and no collision spans the corpora; the window is random, so several budgets
# are tried.  Portable to macOS (no GNU stat).
#
#   ./rehearse_merge.sh
#   CURVE=41 WORKERS=4 DPW=18 LAUNCHLIST="3 4 5 2 6" ./rehearse_merge.sh

set -euo pipefail
cd "$(dirname "$0")"

CURVE=${CURVE:-41}
WORKERS=${WORKERS:-16}
DPW=${DPW:-18}
LAUNCHES=${LAUNCHES:-3}
STEPLIST=${STEPLIST:-"3 2 4"}
ATTEMPTS=${ATTEMPTS:-6}
CLIENT=${CLIENT:-../ecc2k130-cpu}
DIR=${DIR:-$(mktemp -d)}

[ -x "$CLIENT" ] || { echo "build the host client first: (cd .. && make cpu)" >&2; exit 1; }
echo "rehearsal: curve $CURVE, $WORKERS workers, cutoff $DPW, [$STEPLIST] steps x $LAUNCHES launches, $ATTEMPTS attempts, in $DIR"

# Seeds are a function of the run id, so repeating a budget repeats its
# outcome; every attempt uses a fresh run-id range.  With W workers, W-1 of
# every W collisions span two corpora, which is what makes the window usable.
attempt=0
while [ "$attempt" -lt "$ATTEMPTS" ]; do
    for STEPS in $STEPLIST; do
    attempt=$((attempt + 1))
    [ "$attempt" -le "$ATTEMPTS" ] || break
    K=$LAUNCHES
    base=$((attempt * 64))
    rm -rf "$DIR/dp" "$DIR/work"
    mkdir -p "$DIR/dp"
    solo=0
    for i in $(seq 1 "$WORKERS"); do
        # Each worker is one campaign slot: its own run id, its own corpus.
        out=$("$CLIENT" --curve "$CURVE" --threads 1 --steps "$STEPS" --launches "$K" --dp-weight "$DPW" \
                        --verify 0 --run-id "$((base + i))" --dp-file "$DIR/dp/slot-$i.bin" 2>&1)
        if echo "$out" | grep -q "verified \[k\]P == Q"; then solo=1; break; fi
    done
    if [ "$solo" = "1" ]; then
        echo "  steps=$STEPS: a worker solved alone, too much work per worker"
        continue
    fi
    # Detect only, then solve, so both code paths run.
    python3 merge.py --work "$DIR/work" --local "$DIR/dp" --curve "$CURVE" --dp-weight "$DPW" \
        --client "$CLIENT" --buckets 64 --detect-only > "$DIR/detect.json" 2> "$DIR/detect.log"
    n=$(python3 -c "import json;print(json.load(open('$DIR/detect.json'))['collisions'])")
    pts=$(python3 -c "import json;print(json.load(open('$DIR/detect.json'))['corpus'])")
    if [ "$n" = "0" ]; then
        echo "  steps=$STEPS: $pts points, none solved alone, but no collision spanned the corpora"
        continue
    fi
    # A second pass must be incremental: nothing new to ingest, the collision
    # already recorded, so only the solve step should run now.
    python3 merge.py --work "$DIR/work" --local "$DIR/dp" --curve "$CURVE" --dp-weight "$DPW" \
        --client "$CLIENT" > "$DIR/solve.json" 2> "$DIR/solve.log" || true
    if python3 -c "import json,sys; s=json.load(open('$DIR/work/state.json'))['solved']; sys.exit(0 if s and s['verified'] else 1)"; then
        echo "  steps=$STEPS: $WORKERS workers, $pts points, none solved alone"
        echo "  merge found $n cross-corpus collision(s)"
        python3 -c "import json; s=json.load(open('$DIR/work/state.json'))['solved']; print('  k =', s['k'], '(verified; planted match: %s)' % s['matchesPublished'])"
        # The solved state must short-circuit a further pass.
        python3 merge.py --work "$DIR/work" --local "$DIR/dp" --curve "$CURVE" --client "$CLIENT" >/dev/null 2>&1
        echo "REHEARSAL PASSED: cross-corpus collision recovered through merge.py"
        exit 0
    fi
    echo "  steps=$STEPS: $n collision(s) found but none solved:"; tail -20 "$DIR/solve.log"
    exit 1
    done
done
echo "REHEARSAL INCONCLUSIVE: no attempt hit the window" >&2
exit 3

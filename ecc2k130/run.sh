#!/usr/bin/env bash
#
# Launch an ECC2K-130 (or ECC2K-95) search on Modal and keep it running.
#
#   ./run.sh validate            prove the GPU engine before spending anything
#   ./run.sh prewarm             build the Modal image and GPU client for this shape
#   ./run.sh bench               throughput on the challenge curve
#   ./run.sh search              collect, resuming after every container deadline
#   ./run.sh fanout              the same across several GPUs at once
#   ./run.sh merge               scan the corpus for collisions
#   ./run.sh sync                copy new volume DPs into the campaign bucket
#
# A container has a finite life, so a real search is a loop: each pass runs for
# HOURS, is stopped with SIGTERM so the client checkpoints, and the next pass
# resumes from that checkpoint.  Resume only works when the shape matches -- the
# checkpoint header records curve, run id, worker count and batch -- so every
# pass here is launched from the same variables.  Change one between passes and
# the client will refuse the checkpoint and silently start over.
#
# The first Modal run on a cold image can spend tens of minutes compiling CUDA
# while the CLI must keep heartbeating.  search and fanout call validate first
# (PREWARM=1, the default) so that build finishes before the timed pass loop.
# Set PREWARM=0 to skip when the image is already warm.

set -euo pipefail

CURVE=${CURVE:-97}          # 97 = ECC2K-95 (solved 1998, answer known); 131 = ECC2K-130
HOURS=${HOURS:-4}           # per pass; the loop supplies the total
PASSES=${PASSES:-6}         # 0 = until solved
GPU=${ECC_GPU:-RTX-PRO-6000}
BATCH=${BATCH:-8}
THREADS=${THREADS:-128}
LEAF=${LEAF:-17}
WALKS=${WALKS:-4000000}
RUNID=${RUNID:-1}
DPW=${DPW:--1}           # -1 = size the cutoff from the measured rate
LOADMAX=${LOADMAX:-50000000}  # cap the startup corpus reload; 0 = no limit
COUNT=${COUNT:-4}           # fanout width

export ECC_GPU="$GPU"
# Modal relays the container's stdout; keep Python from buffering it so the
# 60-second progress lines arrive while the run is happening, not at the end.
export PYTHONUNBUFFERED=1
cd "$(dirname "$0")"

cmd=${1:-search}
case "$cmd" in
    validate|prewarm|bench|search|fanout|merge|sync) ;;
    *) sed -n '3,11p' "$0"; exit 1 ;;
esac

command -v modal >/dev/null || { echo "modal CLI not found: pip install -U modal" >&2; exit 1; }

shape="--curve $CURVE --batch $BATCH --threads $THREADS --leaf $LEAF --walks $WALKS"
shape="$shape --dp-weight $DPW --load-max $LOADMAX"
build_shape="--batch $BATCH --threads $THREADS --leaf $LEAF"

prewarm_modal() {
    echo "pre-warming Modal image and GPU build on $GPU ($build_shape)..."
    modal run modal_app.py::validate $build_shape
}

case "$cmd" in

validate)
    # Recovers planted discrete logarithms on the GPU itself, so a broken
    # kernel fails in a minute instead of quietly burning a day of credits.
    prewarm_modal
    ;;

prewarm)
    prewarm_modal
    ;;

bench)
    modal run modal_app.py::bench --batch "$BATCH" --threads "$THREADS" --leaf "$LEAF"
    ;;

search)
    if [ "${PREWARM:-1}" != 0 ]; then
        prewarm_modal
    fi
    echo "curve $CURVE on $GPU: $PASSES passes of ${HOURS}h, run id $RUNID"
    echo "checkpoint /data/ckpt/curve$CURVE-run$RUNID.ck, corpus /data/dp/curve$CURVE-run$RUNID.bin"
    pass=1
    while [ "$PASSES" -eq 0 ] || [ "$pass" -le "$PASSES" ]; do
        echo "=== pass $pass ($(date -u +%H:%M:%SZ)) ==="
        # Each pass picks up the checkpoint the last one wrote on its way out.
        # A pass that solves the logarithm prints "k = ..." and exits; stop then
        # rather than starting another.
        if modal run modal_app.py::search \
                $shape --hours "$HOURS" --run-id "$RUNID" | tee /tmp/ecc-pass.$$; then
            if grep -q '"solved": "' /tmp/ecc-pass.$$ && ! grep -q '"solved": null' /tmp/ecc-pass.$$; then
                echo "solved on pass $pass"
                rm -f /tmp/ecc-pass.$$
                if [ "${SYNC:-1}" != 0 ]; then
                    "$0" sync || echo "sync after solve failed; will retry" >&2
                fi
                exit 0
            fi
        else
            echo "pass $pass failed; the checkpoint survives, retrying" >&2
        fi
        rm -f /tmp/ecc-pass.$$
        if [ "${SYNC:-1}" != 0 ]; then
            "$0" sync || echo "sync after pass $pass failed; will retry" >&2
        fi
        pass=$((pass + 1))
    done
    "$0" sync || true
    modal run modal_app.py::merge --curve "$CURVE"
    ;;

fanout)
    # Independent searchers, each with its own run id so their seed spaces stay
    # disjoint, each loading the others' corpora so a cross-worker collision is
    # caught as it happens rather than in the merge afterwards.
    if [ "${PREWARM:-1}" != 0 ]; then
        prewarm_modal
    fi
    echo "curve $CURVE on $COUNT x $GPU: $PASSES passes of ${HOURS}h"
    pass=1
    while [ "$PASSES" -eq 0 ] || [ "$pass" -le "$PASSES" ]; do
        echo "=== pass $pass ($(date -u +%H:%M:%SZ)) ==="
        modal run modal_app.py::fanout $shape --hours "$HOURS" --count "$COUNT" || \
            echo "pass $pass failed; checkpoints survive, retrying" >&2
        pass=$((pass + 1))
    done
    ;;

merge)
    modal run modal_app.py::merge --curve "$CURVE"
    ;;

sync)
    python3 modal_sync.py --curve "$CURVE" --run-id "$RUNID" ${ECC_BUCKET:+--bucket "$ECC_BUCKET"}
    ;;

esac

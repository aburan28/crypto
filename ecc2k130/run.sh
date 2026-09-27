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
#   ./run.sh sync-loop           keep syncing every Modal run until stopped
#   ./run.sh ensure-sync         start sync-loop in tmux when it is not running
#   ./run.sh ingest              copy s3://bucket/dp/ into Postgres (any host)
#   ./run.sh next-run-id         the next free campaign run id on the volume
#   ./run.sh deploy              deploy the app so its functions outlive this shell
#   ./run.sh fleet               keep COUNT campaign runs walking on the deployed app
#
# `fleet` is how the campaign should run. `search`/`fanout` own an ephemeral
# app that Modal stops when this shell's connection drops (and --detach keeps
# only the last spawned function); `fleet` spawns runSearch on the deployed
# app, where calls run to completion on their own, and modal_campaign.py
# re-attaches to them on restart instead of doubling up on a checkpoint.
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
#
# CURVE=131 selects the audited RTX PRO 6000 packed preset (campaign.json):
# batch 16, 256 threads, 385024 workers, dp weight 32, ~14 B iterations/s.
# On that curve RUNID is required and must come from the Modal range
# (8000-9999): a run id an AWS slot has used walks that slot's seeds again,
# and a cutoff other than the campaign's cannot collide with it.  Both are
# refused by modal_app.py and by modal_sync.py; `./run.sh next-run-id` picks
# one, and fanout takes RUNID as the base of COUNT consecutive ids.

set -euo pipefail

CURVE=${CURVE:-97}          # 97 = ECC2K-95 (solved 1998, answer known); 131 = ECC2K-130
HOURS=${HOURS:-4}           # per pass; the loop supplies the total
PASSES=${PASSES:-6}         # 0 = until solved
GPU=${ECC_GPU:-RTX-PRO-6000}
RUNID=${RUNID:-}            # required on curve 131 (see above); 1 elsewhere
COUNT=${COUNT:-4}           # fanout width
CHECKPOINT_EVERY=${CHECKPOINT_EVERY:-60}  # full walk-state file; status uses a 15s header
SYNC_INTERVAL=${SYNC_INTERVAL:-15}        # volume → S3, so the page sees the header
VERIFY=${VERIFY:-0}                       # 0: do not CPU-replay reports (stalls dp32)
# Host threads walking beside each GPU (run id + 1000). Modal bills the cores
# (see modal_app.py CPU_THREADS: ~18 M it/s per AVX-512 thread against
# ~14.8 B for the 6000, about 10x worse per dollar), so this is opt-in.
CPU_THREADS=${CPU_THREADS:-0}
export ECC_CPU_THREADS="$CPU_THREADS"     # sizes the container's CPU request at import

if [ "$CURVE" = 131 ]; then
    # Audited fleet preset — see RTX-PRO6000.md and aws/campaign.json.
    PACKED=${PACKED:-1}
    BATCH=${BATCH:-16}
    THREADS=${THREADS:-256}
    LEAF=${LEAF:-0}
    WORKERS=${WORKERS:-385024}
    WALKS=${WALKS:-$((WORKERS * BATCH))}
    DPW=${DPW:-32}
    LOADMAX=${LOADMAX:-2000000}
    export ECC_CUDA_VERSION="${ECC_CUDA_VERSION:-13.3.1}"
    export ECC_PACKED_SINGLE_PRODUCT="${ECC_PACKED_SINGLE_PRODUCT:-1}"
    export ECC_PACKED_CACHE_DENOM="${ECC_PACKED_CACHE_DENOM:-1}"
    export ECC_PACKED_BY_VALUE="${ECC_PACKED_BY_VALUE:-1}"
    export ECC_PACKED_PERM_SIGMA="${ECC_PACKED_PERM_SIGMA:-3}"
    export ECC_PACKED_POLY_CHAIN="${ECC_PACKED_POLY_CHAIN:-1}"
    export ECC_PACKED_UNROLL_INV="${ECC_PACKED_UNROLL_INV:-1}"
    export ECC_PACKED_PAIR_PRODUCTS="${ECC_PACKED_PAIR_PRODUCTS:-1}"
    export ECC_PACKED_POLY_STATE="${ECC_PACKED_POLY_STATE:-1}"
    export ECC_PACKED_DIRECT_REDUCE="${ECC_PACKED_DIRECT_REDUCE:-1}"
    export ECC_PACKED_GENERATED_PRODUCT="${ECC_PACKED_GENERATED_PRODUCT:-1}"
    export ECC_PACKED_CLMAD="${ECC_PACKED_CLMAD:-1}"
    export ECC_PACKED_STATE_TILE="${ECC_PACKED_STATE_TILE:-256}"
    export ECC_PACKED_WEIGHTED_PREFIX="${ECC_PACKED_WEIGHTED_PREFIX:-2}"
    export ECC_PACKED_COMPACT_STATE="${ECC_PACKED_COMPACT_STATE:-1}"
    export ECC_PACKED_SHARED_SIGMA="${ECC_PACKED_SHARED_SIGMA:-1}"
    export ECC_PACKED_TOP_CLMAD="${ECC_PACKED_TOP_CLMAD:-0}"
    export ECC_PACKED_CLMAD_SQUARE="${ECC_PACKED_CLMAD_SQUARE:-0}"
    export ECC_PACKED_KARAT3="${ECC_PACKED_KARAT3:-0}"
    export ECC_WALK_TABLE="${ECC_WALK_TABLE:-0}"
    export ECC_TABLE_PIVOT_BYTES="${ECC_TABLE_PIVOT_BYTES:-0}"
else
    PACKED=${PACKED:-0}
    BATCH=${BATCH:-8}
    THREADS=${THREADS:-128}
    LEAF=${LEAF:-17}
    WALKS=${WALKS:-4000000}
    DPW=${DPW:--1}
    LOADMAX=${LOADMAX:-50000000}
fi

export ECC_GPU="$GPU"
# Modal relays the container's stdout; keep Python from buffering it so the
# 60-second progress lines arrive while the run is happening, not at the end.
export PYTHONUNBUFFERED=1
cd "$(dirname "$0")"

cmd=${1:-search}
case "$cmd" in
    validate|prewarm|bench|search|fanout|merge|sync|sync-loop|ensure-sync|ingest|next-run-id|deploy|fleet) ;;
    *) sed -n '3,17p' "$0"; exit 1 ;;
esac

command -v modal >/dev/null || { echo "modal CLI not found: pip install -U modal" >&2; exit 1; }

# A pass is launched detached and bounded from this side. Without --detach a
# dropped local connection stops the app: on 2026-09-21 four campaign runners
# were terminated at 15:23Z that way while the local client kept drawing
# "Running (4/4 containers active)" for four hours. The containers stop
# themselves at HOURS; the timeout here is for the client that has stopped
# hearing from them, so the loop can go round and resume the checkpoints.
PASS_TIMEOUT=$(awk "BEGIN{print int($HOURS*3600+1800)}")
modal_pass() {
    timeout --signal=INT --kill-after=60 "$PASS_TIMEOUT" modal run --detach "$@"
}

# The campaign curve refuses an implicit run id: every pass of the loop below
# must launch the same one, so it has to be chosen once, here, by the operator.
require_run_id() {
    if [ "$CURVE" = 131 ]; then
        if [ -z "$RUNID" ]; then
            echo "CURVE=131 needs RUNID in 8000-8999 (./run.sh next-run-id suggests one)" >&2
            exit 1
        fi
        if [ "$RUNID" -lt 8000 ] || [ "$RUNID" -gt 9999 ]; then
            echo "RUNID=$RUNID is outside the Modal campaign range 8000-9999: below it the id" >&2
            echo "belongs to AWS slot $((RUNID - 1)) and the run re-walks that slot's seeds" >&2
            exit 1
        fi
        if [ "$CPU_THREADS" -gt 0 ] && [ "$RUNID" -gt 8999 ]; then
            echo "RUNID=$RUNID cannot carry a CPU walker: its sidecar id $((RUNID + 1000)) is past 9999" >&2
            exit 1
        fi
    else
        RUNID=${RUNID:-1}
    fi
}

packed_flag=
if [ "$PACKED" = 1 ]; then
    packed_flag=--packed
fi

shape="--curve $CURVE --batch $BATCH --threads $THREADS --leaf $LEAF --walks $WALKS"
shape="$shape --dp-weight $DPW --load-max $LOADMAX --verify $VERIFY --cpu-threads $CPU_THREADS"
[ -n "$packed_flag" ] && shape="$shape $packed_flag"
build_shape="--batch $BATCH --threads $THREADS --leaf $LEAF --min-blocks 2"
[ -n "$packed_flag" ] && build_shape="$build_shape $packed_flag"

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
    modal run modal_app.py::bench $build_shape --steps 1024 --launches 32 --workers "${WORKERS:-0}"
    ;;

next-run-id)
    modal run modal_app.py::next_run_id --curve "$CURVE"
    ;;

deploy)
    # The image is baked from this shell's environment (GPU, CPU threads,
    # packed knobs), exactly as `modal run` would bake it.
    modal deploy modal_app.py
    ;;

fleet)
    require_run_id
    "$0" ensure-sync || true
    if [ "${DEPLOY:-1}" != 0 ]; then
        "$0" deploy
    fi
    echo "curve $CURVE on $COUNT x $GPU: ${HOURS}h passes, run ids $RUNID-$((RUNID + COUNT - 1)), cpu threads $CPU_THREADS"
    exec python3 modal_campaign.py --run-id-base "$RUNID" --count "$COUNT" --hours "$HOURS" \
        --passes "$PASSES" $shape --checkpoint-every "$CHECKPOINT_EVERY"
    ;;

search)
    require_run_id
    "$0" ensure-sync || true
    if [ "${PREWARM:-1}" != 0 ]; then
        prewarm_modal
    fi
    echo "curve $CURVE on $GPU: $PASSES passes of ${HOURS}h, run id $RUNID"
    if [ "$PACKED" = 1 ]; then
        echo "packed preset: batch=$BATCH threads=$THREADS workers=$WORKERS walks=$WALKS dp-weight=$DPW"
    fi
    echo "checkpoint /data/ckpt/curve$CURVE-run$RUNID.ck every ${CHECKPOINT_EVERY}s, corpus /data/dp/curve$CURVE-run$RUNID.bin"
    pass=1
    while [ "$PASSES" -eq 0 ] || [ "$pass" -le "$PASSES" ]; do
        echo "=== pass $pass ($(date -u +%H:%M:%SZ)) ==="
        # Each pass picks up the checkpoint the last one wrote on its way out.
        # A pass that solves the logarithm prints "k = ..." and exits; stop then
        # rather than starting another.
        if modal_pass modal_app.py::search \
                $shape --hours "$HOURS" --run-id "$RUNID" \
                --checkpoint-every "$CHECKPOINT_EVERY" | tee /tmp/ecc-pass.$$; then
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
    # Independent searchers on run ids RUNID .. RUNID+COUNT-1, so their seed
    # spaces stay disjoint, each loading the others' corpora so a cross-worker
    # collision is caught as it happens rather than in the merge afterwards.
    require_run_id
    "$0" ensure-sync || true
    if [ "${PREWARM:-1}" != 0 ]; then
        prewarm_modal
    fi
    echo "curve $CURVE on $COUNT x $GPU: $PASSES passes of ${HOURS}h, run ids $RUNID-$((RUNID + COUNT - 1))"
    pass=1
    while [ "$PASSES" -eq 0 ] || [ "$pass" -le "$PASSES" ]; do
        echo "=== pass $pass ($(date -u +%H:%M:%SZ)) ==="
        modal_pass modal_app.py::fanout $shape --hours "$HOURS" --count "$COUNT" \
            --run-id-base "$RUNID" --checkpoint-every "$CHECKPOINT_EVERY" || \
            echo "pass $pass failed; checkpoints survive, retrying" >&2
        pass=$((pass + 1))
    done
    ;;

merge)
    modal run modal_app.py::merge --curve "$CURVE"
    ;;

sync)
    require_run_id
    python3 modal_sync.py --curve "$CURVE" --run-id "$RUNID" ${ECC_BUCKET:+--bucket "$ECC_BUCKET"}
    ;;

sync-loop)
    # Discover every curve*-run*.bin on the Modal volume each pass so fanout
    # workers need no sync configuration. Failures are logged and retried.
    python3 modal_sync.py --curve "$CURVE" --all-runs --watch "$SYNC_INTERVAL" \
        ${SYNC_RUN_IDS:+--run-ids "$SYNC_RUN_IDS"} \
        ${ECC_BUCKET:+--bucket "$ECC_BUCKET"}
    ;;

ensure-sync)
    SESSION=${MODAL_SYNC_SESSION:-ecc2k130-modal-sync}
    LOG=${MODAL_SYNC_LOG:-/tmp/ecc2k130-modal-sync.log}
    TMUX=(tmux -f /exec-daemon/tmux.portal.conf)
    if "${TMUX[@]}" has-session -t "=$SESSION" 2>/dev/null; then
        echo "sync loop already running in tmux session $SESSION"
        exit 0
    fi
    "${TMUX[@]}" new-session -d -s "$SESSION" -c "$(dirname "$0")" -- "${SHELL:-bash}" -l -c \
        "while true; do ./run.sh sync-loop 2>&1 || sleep 30; done | tee -a $LOG"
    echo "started sync loop in tmux session $SESSION (log: $LOG)"
    ;;

ingest)
    shift
    INGEST_ENSURE_ACCESS=${INGEST_ENSURE_ACCESS:-1} exec aws/ingest.sh "$@"
    ;;

esac

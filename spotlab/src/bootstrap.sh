#!/bin/bash
# spotlab bootstrap: runs once as root at first boot (EC2 user-data, GCE
# startup-script). Fetch the agent, run one attempt, upload the agent log,
# power off. The instance is configured to terminate (EC2) or be deleted
# (GCE) when it powers off, so nothing is left running or billed.
set -u
export SPOTLAB_STORE='@STORE@'
export SPOTLAB_JOB_ID='@JOB@'
export SPOTLAB_ATTEMPT='@ATTEMPT@'
export SPOTLAB_PROVIDER='@PROVIDER@'

# Hard ceiling, set before anything that could hang.
shutdown -h +@DEADLINE_MIN@ 'spotlab: job time budget exhausted' || true

LOG=/var/log/spotlab-agent.log
exec >>"$LOG" 2>&1
echo "spotlab bootstrap: job $SPOTLAB_JOB_ID attempt $SPOTLAB_ATTEMPT on $SPOTLAB_PROVIDER, $(date -u)"

# --- provider-specific setup from the config (extra_bootstrap) ---
@EXTRA@
# ---

fetch() {
  case "$SPOTLAB_STORE" in
    s3://*) aws s3 cp --only-show-errors "$SPOTLAB_STORE/$1" "$2" ;;
    gs://*) gcloud storage cp "$SPOTLAB_STORE/$1" "$2" ;;
    *) echo "unsupported store $SPOTLAB_STORE"; return 1 ;;
  esac
}
put() {
  case "$SPOTLAB_STORE" in
    s3://*) aws s3 cp --only-show-errors "$1" "$SPOTLAB_STORE/$2" ;;
    gs://*) gcloud storage cp "$1" "$SPOTLAB_STORE/$2" ;;
  esac
}

dir=/opt/spotlab
mkdir -p "$dir"
ok=0
for i in 1 2 3 4 5; do
  if fetch "bin/spotlab-$(uname -m)" "$dir/spotlab"; then ok=1; break; fi
  sleep $((i * 3))
done
if [ "$ok" = 1 ]; then
  chmod +x "$dir/spotlab"
  "$dir/spotlab" agent --store "$SPOTLAB_STORE" --job "$SPOTLAB_JOB_ID" \
    --attempt "$SPOTLAB_ATTEMPT" --provider "$SPOTLAB_PROVIDER" \
    --workdir "/var/lib/spotlab/$SPOTLAB_JOB_ID/a$SPOTLAB_ATTEMPT"
  echo "spotlab bootstrap: agent exited $?"
else
  echo "spotlab bootstrap: could not fetch the agent binary"
fi
put "$LOG" "jobs/$SPOTLAB_JOB_ID/attempts/$(printf %04d "$SPOTLAB_ATTEMPT")/agent.log" || true
shutdown -h now

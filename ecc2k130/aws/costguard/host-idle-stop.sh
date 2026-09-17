#!/usr/bin/env bash
#
# Optional host-local idle/age stop for durable GPU boxes.
# Defaults are OFF (0) — the account monthly budget in costguard.sh is the
# primary runaway control. Enable only if you want a forgotten Cursor box to
# stop itself without waiting for the monthly cap:
#
#   COSTGUARD_MAX_AGE_HOURS=168 COSTGUARD_IDLE_HOURS=48 sudo -E ./install-host.sh
#
# Age is measured from when costguard was armed (not boot). Keepalive refreshes
# the age clock. Touch /var/lib/costguard/exempt to disable.

set -euo pipefail

STATE_DIR=${COSTGUARD_STATE_DIR:-/var/lib/costguard}
LOG=${COSTGUARD_LOG:-/var/log/costguard-idle-stop.log}
MAX_AGE_HOURS=${COSTGUARD_MAX_AGE_HOURS:-0}
IDLE_HOURS=${COSTGUARD_IDLE_HOURS:-0}
IDLE_LOAD=${COSTGUARD_IDLE_LOAD:-0.3}
READYZ_URL=${COSTGUARD_READYZ_URL:-http://127.0.0.1:9182/readyz}
IDLE_MARKER=$STATE_DIR/idle-since
KEEPALIVE=$STATE_DIR/keepalive
EXEMPT=$STATE_DIR/exempt
ARMED_AT=$STATE_DIR/armed-at

mkdir -p "$STATE_DIR"
chmod 755 "$STATE_DIR"

log() { printf '%s %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$LOG" >/dev/null; }

if [[ -f "$EXEMPT" ]]; then
  log "exempt marker present; no action"
  exit 0
fi

if [[ ! -f "$ARMED_AT" ]]; then
  date +%s > "$ARMED_AT"
  touch "$KEEPALIVE"
  log "armed costguard; created armed-at and keepalive"
fi

now=$(date +%s)
armed_at=$(cat "$ARMED_AT")
age_hours=$(python3 -c "print(max(0, ($now - int('$armed_at')) / 3600.0))")
load1=$(awk '{print $1}' /proc/loadavg)

if python3 -c "import sys; sys.exit(0 if float('$MAX_AGE_HOURS') > 0 else 1)"; then
  keepalive_ok=0
  if [[ -f "$KEEPALIVE" ]]; then
    age_sec=$(( now - $(stat -c %Y "$KEEPALIVE") ))
    if (( age_sec < MAX_AGE_HOURS * 3600 )); then
      keepalive_ok=1
    fi
  fi
  if python3 - "$age_hours" "$MAX_AGE_HOURS" "$keepalive_ok" <<'PY'
import sys
age, limit, keepalive = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
sys.exit(0 if (age >= limit and keepalive == 0) else 1)
PY
  then
    log "max age reached: armed_age_h=$age_hours limit=$MAX_AGE_HOURS keepalive=$keepalive_ok; shutting down"
    /sbin/shutdown -h now "costguard: max age ${MAX_AGE_HOURS}h"
    exit 0
  fi
fi

if ! python3 -c "import sys; sys.exit(0 if float('$IDLE_HOURS') > 0 else 1)"; then
  log "idle stop disabled (COSTGUARD_IDLE_HOURS=0); age_h=$age_hours load=$load1"
  exit 0
fi

claimed=unknown
if out=$(curl -fsS --max-time 2 "$READYZ_URL" 2>/dev/null); then
  if echo "$out" | grep -q '"claimed":true'; then
    claimed=true
  elif echo "$out" | grep -q '"claimed":false'; then
    claimed=false
  fi
fi

busy=0
if python3 - "$load1" "$IDLE_LOAD" <<'PY'
import sys
sys.exit(0 if float(sys.argv[1]) >= float(sys.argv[2]) else 1)
PY
then
  busy=1
fi
[[ "$claimed" == "true" ]] && busy=1

if [[ "$busy" -eq 1 ]]; then
  rm -f "$IDLE_MARKER"
  log "busy load=$load1 claimed=$claimed; idle timer cleared"
  exit 0
fi

if [[ ! -f "$IDLE_MARKER" ]]; then
  echo "$now" > "$IDLE_MARKER"
  log "became idle load=$load1 claimed=$claimed; starting idle timer"
  exit 0
fi

idle_since=$(cat "$IDLE_MARKER")
idle_hours=$(python3 -c "print(($now - int('$idle_since')) / 3600.0)")
if python3 - "$idle_hours" "$IDLE_HOURS" <<'PY'
import sys
sys.exit(0 if float(sys.argv[1]) >= float(sys.argv[2]) else 1)
PY
then
  log "idle for ${idle_hours}h (limit ${IDLE_HOURS}h) load=$load1 claimed=$claimed; shutting down"
  /sbin/shutdown -h now "costguard: idle ${IDLE_HOURS}h"
  exit 0
fi

log "idle for ${idle_hours}h < ${IDLE_HOURS}h; waiting"

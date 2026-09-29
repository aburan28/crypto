#!/bin/zsh
# LaunchAgent entry for the ECC2K-130 distinguished-point ingester.
# Installed by install-macos.sh; launchd does not inherit a login shell.
set -eu
umask 077

if [[ -z ${HOME:-} ]]; then
  print -u2 "HOME is unset; this runner is a per-user LaunchAgent"
  exit 1
fi

export PATH=/Library/Frameworks/Python.framework/Versions/3.13/bin:/opt/homebrew/bin:/usr/bin:/bin:/usr/sbin:/sbin
export AWS_PROFILE=${DP_INGEST_AWS_PROFILE:-ecc2k130}
export AWS_DEFAULT_REGION=${AWS_DEFAULT_REGION:-us-west-2}
unset AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY AWS_SESSION_TOKEN

ingest_root="$HOME/Library/Application Support/ECC2K130/ingest"
mkdir -p "$ingest_root"
print -r -- "$$" > "$ingest_root/launchd.pid"
print -r -- "$(date -u +%Y-%m-%dT%H:%M:%SZ) launchd starting dp_ingest"

# ingest.sh opens this host's /32 only when it starts. A laptop egress
# address can change while that process stays up, and every later connection
# to rho-dp then times out. Refresh the rule for the life of the job.
refresh_access() {
  /bin/zsh "$ingest_root/ingest.sh" ensure-access
}

refresh_access || print -u2 "ensure-access failed at start; ingest will retry"

(
  while true; do
    sleep 120
    refresh_access || print -u2 "ensure-access failed; will retry"
  done
) &
refresher=$!

cleanup() {
  kill $refresher 2>/dev/null || true
  if [[ -n ${child:-} ]]; then
    kill $child 2>/dev/null || true
  fi
}
trap cleanup EXIT INT TERM

INGEST_ENSURE_ACCESS=0 /bin/zsh "$ingest_root/ingest.sh" &
child=$!
rc=0
wait $child || rc=$?
exit $rc

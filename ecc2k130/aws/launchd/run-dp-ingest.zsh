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

# A laptop's egress address changes between launches. ingest.sh adds this
# host's /32 only when this is set, and only at process start.
export INGEST_ENSURE_ACCESS=1

ingest_root="$HOME/Library/Application Support/ECC2K130/ingest"
mkdir -p "$ingest_root"
print -r -- "$$" > "$ingest_root/launchd.pid"
print -r -- "$(date -u +%Y-%m-%dT%H:%M:%SZ) launchd starting dp_ingest"

exec /bin/zsh "$ingest_root/ingest.sh"

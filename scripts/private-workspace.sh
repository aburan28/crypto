#!/usr/bin/env bash
# Preinstalled-only startup for an isolated Conductor/Cairn workspace.
# No package installation, release discovery, seed refresh, or tunnel mutation.
set -euo pipefail

fail() { echo "private-workspace: $*" >&2; exit 1; }
action=${1:-check}
case "$action" in check|run) ;; *) fail "usage: $0 [check|run]" ;; esac

for binary in conductor conductor-mcp cairn; do
  command -v "$binary" >/dev/null || fail "preinstall $binary from your approved offline bundle"
done

# A legacy binary silently ignores unknown environment variables. Require the
# native profile's help marker before permitting it to open any network sockets.
env -u CONDUCTOR_TOKEN cairn --version 2>/dev/null | grep embedded >/dev/null \
  || fail "preinstall Cairn with its embedded UI; offline mode never downloads a replacement"
env -u CONDUCTOR_TOKEN cairn run --help 2>&1 | grep -F 'CAIRN_NETWORK=private' >/dev/null \
  || fail "this Cairn binary does not advertise the private-network profile; install a compatible approved release"

[[ -n "${CONDUCTOR_ENDPOINT:-}" && -n "${CONDUCTOR_TOKEN:-}" ]] \
  || fail "configure CONDUCTOR_ENDPOINT and a dedicated CONDUCTOR_TOKEN for the internal control plane"
[[ "$CONDUCTOR_ENDPOINT" =~ ^(https?)://([A-Za-z0-9][A-Za-z0-9.-]*|\[[0-9A-Fa-f:]+\])(:([0-9]{1,5}))?/?$ ]] \
  || fail "Conductor endpoint must be an HTTP(S) origin; put credentials only in CONDUCTOR_TOKEN"
scheme=${BASH_REMATCH[1]}
host=${BASH_REMATCH[2]}
port=${BASH_REMATCH[4]:-}
[[ -z "$port" ]] || ((10#$port >= 1 && 10#$port <= 65535)) || fail "Conductor endpoint has an invalid port"
[[ "$scheme" == https || "$host" == localhost || "$host" == 127.0.0.1 || "$host" == '[::1]' ]] \
  || fail "Conductor requires HTTPS, except for an explicit loopback origin"

bootstrap=${CAIRN_PRIVATE_BOOTSTRAP:-}
[[ -z "$bootstrap" || -f "$bootstrap" && -r "$bootstrap" ]] \
  || fail "CAIRN_PRIVATE_BOOTSTRAP must be a readable local bootstrap file"
[[ -z "${CAIRN_PRIVATE_ROOT:-}" || -d "$CAIRN_PRIVATE_ROOT" ]] \
  || fail "CAIRN_PRIVATE_ROOT must name the local, preinstalled verifier bundle"

export CAIRN_NETWORK=private
export CAIRN_SEEDS=off CAIRN_PEERS=bootstrap
export CAIRN_BEACON_PORT=off CAIRN_PORTMAP=off CAIRN_PORTMAP_HTTP=off

if [[ "$action" == check ]]; then
  echo "private-workspace: compatible local tools found; automatic installs and public discovery are disabled"
  exit 0
fi

state_dir=${CAIRN_PRIVATE_STATE:-$HOME/.cairn/private-node}
mkdir -p "$state_dir"
chmod 700 "$state_dir"
export CAIRN_DATA="$state_dir"
# A relay does not start an autonomous job executor or spend on objectives.
export CAIRN_ROLES=relay CAIRN_MCP_MAX_SPEND=0
verifier_root=${CAIRN_PRIVATE_ROOT:-$state_dir}
args=(--root "$verifier_root" run --no-mcp
      --listen "${CAIRN_NODE_P2P:-127.0.0.1:9100}"
      --serve "${CAIRN_NODE_HTTP:-127.0.0.1:9180}")
[[ -z "$bootstrap" ]] || args+=(--bootstrap "$bootstrap")
# The relay has its own identity; the team's control-plane credential is not
# authority it needs. This only changes the child, not the caller's shell.
unset CONDUCTOR_TOKEN
exec cairn "${args[@]}"

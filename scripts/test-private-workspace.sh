#!/usr/bin/env bash
# Native shell contract checks; external/network/install commands are tripwires.
set -euo pipefail
repo=$(cd "$(dirname "$0")/.." && pwd)
fixture=$(mktemp -d)
trap 'rm -rf "$fixture"' EXIT
mkdir "$fixture/bin"
export PRIVATE_TEST_TRACE="$fixture/trace"
export PRIVATE_TEST_RUN="$fixture/run"
unset CAIRN_PRIVATE_BOOTSTRAP CAIRN_PRIVATE_ROOT CAIRN_NODE_P2P CAIRN_NODE_HTTP
unset PRIVATE_TEST_HEADLESS PRIVATE_TEST_LEGACY PRIVATE_TEST_VERSION_STATUS PRIVATE_TEST_HELP_STATUS PRIVATE_TEST_RUN_STATUS
export HOME="$fixture/home"
mkdir "$HOME"
export PATH="$fixture/bin:$PATH"
export CONDUCTOR_ENDPOINT=https://conductor.internal:8443
export CONDUCTOR_TOKEN=test-only-not-a-real-credential
export WORKSPACE_OFFLINE=1
export CRYPTO_REPO_DIR="$repo" CLAUDE_PROJECT_DIR="$repo"
export CAIRN_PRIVATE_STATE="$fixture/state"

cat > "$fixture/bin/cairn" <<'STUB'
#!/usr/bin/env bash
set -euo pipefail
[[ -z "${CONDUCTOR_TOKEN+x}" ]] || { echo 'Conductor credential leaked into Cairn' >&2; exit 1; }
case "${1:-}" in
  --version)
    echo 'cairn fixture'
    if [[ "${PRIVATE_TEST_HEADLESS:-0}" == 1 ]]; then
      echo '  ui       not embedded -- cairn run refuses'
    else
      echo '  ui       embedded -- cairn run will start'
    fi
    exit "${PRIVATE_TEST_VERSION_STATUS:-0}" ;;
  run)
    if [[ "${2:-}" == --help ]]; then
      [[ "${PRIVATE_TEST_LEGACY:-0}" == 1 ]] || echo 'CAIRN_NETWORK=private'
      exit "${PRIVATE_TEST_HELP_STATUS:-0}"
    fi ;;
esac
[[ "$CAIRN_NETWORK" == private && "$CAIRN_SEEDS" == off && "$CAIRN_PEERS" == bootstrap ]]
[[ "$CAIRN_PORTMAP" == off && "$CAIRN_PORTMAP_HTTP" == off && "$CAIRN_BEACON_PORT" == off ]]
[[ "$CAIRN_ROLES" == relay && "$CAIRN_MCP_MAX_SPEND" == 0 ]]
[[ "$CAIRN_DATA" == "$CAIRN_PRIVATE_STATE" ]]
printf '%s\n' "$@" > "$PRIVATE_TEST_RUN"
exit "${PRIVATE_TEST_RUN_STATUS:-0}"
STUB
cat > "$fixture/bin/conductor" <<'STUB'
#!/usr/bin/env bash
if [[ "${1:-}" == --version ]]; then
  echo 'conductor fixture'
  exit 0
fi
echo "${0##*/} service/auth invocation" >> "$PRIVATE_TEST_TRACE"
exit 90
STUB
cp "$fixture/bin/conductor" "$fixture/bin/conductor-mcp"
cat > "$fixture/bin/tripwire" <<'STUB'
#!/usr/bin/env bash
echo "${0##*/}" >> "$PRIVATE_TEST_TRACE"
exit 90
STUB
chmod +x "$fixture/bin/"*
for binary in curl wget git go cargo npm uv pip pip3 apt-get apt sudo brew docker runuser postgres pg_ctl initdb psql systemctl service tailscale ssh; do
  ln -s tripwire "$fixture/bin/$binary"
done

cases=0
expect_failure() {
  cases=$((cases + 1))
  if "$@" > "$fixture/out" 2> "$fixture/err"; then
    echo "Unexpected success: $*" >&2; exit 1
  fi
}

bash "$repo/scripts/private-workspace.sh" check > "$fixture/out"
bash "$repo/.claude/hooks/session-start.sh" > "$fixture/out"
bash "$repo/scripts/setup-cloud-environment.sh" > "$fixture/out"
bash "$repo/scripts/setup-workstation.sh" > "$fixture/out"
[[ ! -e "$CAIRN_PRIVATE_STATE" && ! -e "$PRIVATE_TEST_RUN" ]]

expect_failure env PRIVATE_TEST_LEGACY=1 bash "$repo/scripts/private-workspace.sh" check
expect_failure env PRIVATE_TEST_LEGACY=1 bash "$repo/scripts/setup-cloud-environment.sh"
expect_failure env PRIVATE_TEST_HEADLESS=1 bash "$repo/scripts/private-workspace.sh" check
expect_failure env PRIVATE_TEST_HEADLESS=1 bash "$repo/scripts/setup-workstation.sh"
expect_failure env CONDUCTOR_TOKEN= bash "$repo/.claude/hooks/session-start.sh"
expect_failure env CONDUCTOR_ENDPOINT=http://10.77.0.1:8080 bash "$repo/scripts/private-workspace.sh" check
expect_failure env CONDUCTOR_ENDPOINT=https://user:secret@controller.internal bash "$repo/scripts/private-workspace.sh" check
expect_failure env CONDUCTOR_ENDPOINT=https:// bash "$repo/scripts/private-workspace.sh" check
expect_failure env CONDUCTOR_ENDPOINT=https://conductor.internal:99999 bash "$repo/scripts/private-workspace.sh" check
expect_failure env CONDUCTOR_ENDPOINT=https://conductor.internal/api bash "$repo/scripts/private-workspace.sh" check
env CONDUCTOR_ENDPOINT=http://127.0.0.1:8080 bash "$repo/scripts/private-workspace.sh" check > "$fixture/out"
env 'CONDUCTOR_ENDPOINT=https://[fd00::1]:8443' bash "$repo/scripts/private-workspace.sh" check > "$fixture/out"
mv "$fixture/bin/conductor-mcp" "$fixture/bin/missing-mcp"
expect_failure env PATH="$fixture/bin" /bin/bash "$repo/scripts/private-workspace.sh" check
mv "$fixture/bin/missing-mcp" "$fixture/bin/conductor-mcp"
expect_failure env CAIRN_PRIVATE_BOOTSTRAP="$fixture/missing.json" bash "$repo/scripts/private-workspace.sh" check
for entry in .claude/hooks/session-start.sh scripts/setup-cloud-environment.sh scripts/setup-workstation.sh; do
  expect_failure env WORKSPACE_OFFLINE=typo bash "$repo/$entry"
done

# Failed probes must not count as compatible even when they print the marker.
expect_failure env PRIVATE_TEST_VERSION_STATUS=17 bash "$repo/scripts/private-workspace.sh" check
expect_failure env PRIVATE_TEST_HELP_STATUS=18 bash "$repo/scripts/private-workspace.sh" check
for binary in conductor conductor-mcp cairn; do
  mv "$fixture/bin/$binary" "$fixture/bin/hidden-$binary"
  expect_failure env PATH="$fixture/bin" /bin/bash "$repo/scripts/private-workspace.sh" check
  grep -Fq "preinstall $binary" "$fixture/err"
  mv "$fixture/bin/hidden-$binary" "$fixture/bin/$binary"
done
for endpoint in 'http://localhost.evil' 'http://127.0.0.2' 'https://host:0' 'https://host:00000' 'https://host:65536' 'https://host?token=secret' 'https://host/#fragment' 'ftp://host' 'https://host\path'; do
  expect_failure env CONDUCTOR_ENDPOINT="$endpoint" bash "$repo/scripts/private-workspace.sh" check
done
for endpoint in 'http://localhost:1/' 'http://[::1]:65535' 'https://controller.internal:65535/'; do
  env CONDUCTOR_ENDPOINT="$endpoint" bash "$repo/scripts/private-workspace.sh" check > "$fixture/out"
done
expect_failure env CAIRN_PRIVATE_ROOT="$fixture/missing-root" bash "$repo/scripts/private-workspace.sh" run
expect_failure env CAIRN_PRIVATE_BOOTSTRAP="$fixture" bash "$repo/scripts/private-workspace.sh" run
expect_failure bash "$repo/scripts/private-workspace.sh" invalid-action
grep -Fq 'usage:' "$fixture/err"
[[ ! -e "$CAIRN_PRIVATE_STATE" && ! -e "$PRIVATE_TEST_RUN" ]]

printf '{}\n' > "$fixture/internal bootstrap.json"
export CAIRN_PRIVATE_BOOTSTRAP="$fixture/internal bootstrap.json"
export CAIRN_NODE_P2P=10.77.0.2:9100 CAIRN_NODE_HTTP=10.77.0.2:9180
# Ambient public settings must not weaken the native private launch.
export CAIRN_SEEDS=builtin CAIRN_PEERS=open CAIRN_PORTMAP=on CAIRN_PORTMAP_HTTP=on
bash "$repo/scripts/private-workspace.sh" run
cat > "$fixture/expected" <<ARGS
--root
$CAIRN_PRIVATE_STATE
run
--no-mcp
--listen
$CAIRN_NODE_P2P
--serve
$CAIRN_NODE_HTTP
--bootstrap
$CAIRN_PRIVATE_BOOTSTRAP
ARGS
cmp "$fixture/expected" "$PRIVATE_TEST_RUN"
[[ "$CONDUCTOR_TOKEN" == test-only-not-a-real-credential ]]
# Reused state is protected and a verifier path with spaces remains one arg.
chmod 755 "$CAIRN_PRIVATE_STATE"
mkdir "$fixture/verifier bundle"
env CAIRN_PRIVATE_ROOT="$fixture/verifier bundle" bash "$repo/scripts/private-workspace.sh" run
[[ $(sed -n '2p' "$PRIVATE_TEST_RUN") == "$fixture/verifier bundle" ]]
[[ $(stat -c %a "$CAIRN_PRIVATE_STATE") == 700 ]]
unset CAIRN_PRIVATE_BOOTSTRAP CAIRN_NODE_P2P CAIRN_NODE_HTTP
bash "$repo/scripts/private-workspace.sh" run
cat > "$fixture/expected" <<ARGS
--root
$CAIRN_PRIVATE_STATE
run
--no-mcp
--listen
127.0.0.1:9100
--serve
127.0.0.1:9180
ARGS
cmp "$fixture/expected" "$PRIVATE_TEST_RUN"
# Native failure is propagated; no installer/retry/fallback hides the failure.
set +e
env PRIVATE_TEST_RUN_STATUS=37 bash "$repo/scripts/private-workspace.sh" run
status=$?
set -e
[[ "$status" == 37 ]]
[[ $(stat -c %a "$CAIRN_PRIVATE_STATE") == 700 ]]
[[ ! -e "$PRIVATE_TEST_TRACE" ]] || { cat "$PRIVATE_TEST_TRACE" >&2; exit 1; }
echo 'private workspace checks passed: preinstalled-only entrypoints, fail-closed legacy/config handling, private launch, zero install/network tripwires'
echo "$cases rejection cases passed"

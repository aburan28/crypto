#!/bin/bash
# SessionStart hook for Claude Code on the web. It installs Conductor and cairn
# and brings up a Conductor control plane, so the Conductor integration #1422
# committed (.mcp.json, the pre-edit hooks, the CLAUDE.md block) and the cairn
# MCP server in .mcp.json work in every cloud session. docs/CONDUCTOR.md
# explains each step.
#
# It starts no autonomous agent: no `conductor worker`, no `cairn agent run`.
#
# Idempotent: binaries are installed only when missing, so a container restored
# from the cached image goes straight to starting the control plane.
set -euo pipefail

if [ "${CLAUDE_CODE_REMOTE:-}" != "true" ]; then
  exit 0
fi

state="$HOME/.local/state/workspace-tools"
mkdir -p "$state" "$HOME/.cairn"
chmod 700 "$state"
# The hook's stdout becomes session context, so keep it for one summary line and
# send everything else to a private log: `conductor up` prints a login token.
exec 3>&1
exec >>"$state/session-start.log" 2>&1
chmod 600 "$state/session-start.log"
echo "== $(date -u +%FT%TZ) session start"

project=crypto
repo="${CLAUDE_PROJECT_DIR:-$(cd "$(dirname "$0")/../.." && pwd)}"

# Conductor, through the Go module proxy. The release installer cannot download
# here (the session's GitHub proxy answers 403), and @latest resolves to v0.1.0,
# whose go.mod declares the wrong module path.
if ! command -v conductor >/dev/null || ! command -v conductor-mcp >/dev/null; then
  GOBIN=/usr/local/bin GOFLAGS=-modcacherw \
    go install "github.com/aburan28/conductor/cmd/...@${CONDUCTOR_VERSION:-main}"
fi

# cairn, built from its public repository for the same reason. The `ui` feature
# stays off: `cairn mcp` does not need it.
if ! command -v cairn >/dev/null; then
  cargo install --locked --root /usr/local \
    --git https://github.com/aburan28/cairn ${CAIRN_REV:+--rev "$CAIRN_REV"} cairn
fi

# A shared control plane when the environment names one; otherwise one in this
# container on the image's Postgres 16, which coordinates only this container.
if [ -n "${CONDUCTOR_ENDPOINT:-}" ] && [ -n "${CONDUCTOR_TOKEN:-}" ]; then
  plane="the shared control plane at $CONDUCTOR_ENDPOINT"
else
  pgbin=/usr/lib/postgresql/16/bin
  pgdata=/var/lib/conductor-pg
  port=55432
  if [ ! -s "$pgdata/PG_VERSION" ]; then
    install -d -o postgres -g postgres -m 700 "$pgdata"
    runuser -u postgres -- "$pgbin/initdb" -D "$pgdata" -U conductor --auth=trust -E UTF8
  fi
  if ! "$pgbin/pg_isready" -q -h 127.0.0.1 -p "$port"; then
    # A container restored from the cached image keeps the postmaster.pid of a
    # server that no longer runs, and Postgres refuses to start while that pid
    # names any live process. Remove it unless the pid is a live postgres.
    if [ -f "$pgdata/postmaster.pid" ]; then
      pid="$(head -n 1 "$pgdata/postmaster.pid")"
      if [ "$(cat "/proc/$pid/comm" 2>/dev/null)" != postgres ] \
        || grep -q '^State:[[:space:]]*Z' "/proc/$pid/status" 2>/dev/null; then
        rm -f "$pgdata/postmaster.pid"
      fi
    fi
    runuser -u postgres -- "$pgbin/pg_ctl" -D "$pgdata" -l "$pgdata/server.log" -w \
      -o "-p $port -k /tmp -c listen_addresses=127.0.0.1" start
  fi
  if ! "$pgbin/psql" -h 127.0.0.1 -p "$port" -U conductor -d postgres -tAc \
      "select 1 from pg_database where datname = 'conductor'" | grep -q 1; then
    "$pgbin/createdb" -h 127.0.0.1 -p "$port" -U conductor conductor
  fi
  (cd "$repo" && conductor up --project "$project" \
    --dsn "postgres://conductor@127.0.0.1:$port/conductor?sslmode=disable")
  plane="a control plane in this container (http://localhost:8080; it coordinates nothing outside it)"
fi

conductor_version="$(conductor version 2>/dev/null | head -n 1)"
cairn_version="$(cairn --version 2>/dev/null | head -n 1)"
echo "ready: $conductor_version; $cairn_version; $plane"
echo "Workspace tools: $conductor_version and $cairn_version are installed; Conductor is using $plane. Log: $state/session-start.log" >&3

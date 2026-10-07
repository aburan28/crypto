# Conductor

This repository is wired to [Conductor](https://github.com/aburan28/conductor),
the coordination control plane that tells concurrent sessions, human or agent,
who holds which files. #1422 committed what `conductor integrate claude
--project crypto` generates:

| File | What it does |
| --- | --- |
| `.mcp.json` | starts the `conductor` MCP server (`conductor-mcp --project crypto`) |
| `.claude/settings.json` | `PreToolUse` checks every Edit/Write and reserves the path; `SessionStart`/`SessionEnd` register the session; `Stop`/`PreCompact` save a checkpoint |
| `CLAUDE.md` | the "check before you edit" block between the `conductor:begin`/`end` markers |

None of it does anything until the three binaries (`conductor`,
`conductor-mcp`, `conductord`) are on `PATH` and a control plane is reachable.
Without them the hooks print one warning and let the edit through, the MCP
server fails to start (`ENOENT: conductor-mcp`), and the `conductor check`
that `CLAUDE.md` asks for cannot run. This page is how to install them.

## On your own machine

You need git and curl, plus PostgreSQL 16+ or Docker for the control plane.
Go 1.25.14+ is needed only to build from source; an older `go` downloads the
right toolchain itself. Pick one way to install:

```sh
# A release build (macOS/Linux, amd64/arm64), checked against SHA256SUMS,
# into ~/.local/bin. Put that directory on PATH.
curl -fsSL https://raw.githubusercontent.com/aburan28/conductor/main/scripts/install-release.sh \
  | bash -s -- aburan28/conductor ~/.local/bin

# From a clone, building that checkout; adds ~/.local/bin to your shell's PATH.
git clone https://github.com/aburan28/conductor && cd conductor && make install-local

# With Go. Use @main: @latest resolves to the v0.1.0 tag, whose go.mod still
# declares github.com/adamburan/conductor, and `go install` refuses it.
go install github.com/aburan28/conductor/cmd/...@main
```

Then, from the root of this repository:

```sh
conductor up       # Postgres (Docker only if none answers), the control plane on 127.0.0.1:8080, your login
conductor doctor   # control plane, database, binaries, and which coding tools are connected
```

Do not run `conductor integrate claude` again: its output is already
committed. To see what a newer Conductor would change, use
`conductor integrate claude --print`, and commit any change the same way #1422
did.

`.conductor/`, which holds the project policy and the `WORKFLOW.md` that
`CLAUDE.md` tells agents to read, is not in the repository yet.
`conductor init` scaffolds it. It also adds the same managed block to
`AGENTS.md`, so whether to commit that is a decision for the repository
owner, not a side effect of installing.

Start sessions with `conductor wrap claude`, or run `claude` directly and let
the hooks register it.

## In Claude Code on the web

A cloud session runs in a fresh container, so nothing is installed until
something installs it. The repository does that itself:
`.claude/hooks/session-start.sh` is registered as a `SessionStart` hook in
`.claude/settings.json`, and in a cloud session (`CLAUDE_CODE_REMOTE=true`)
it does three things.

- **Installs Conductor.** It runs
  `GOBIN=/usr/local/bin go install github.com/aburan28/conductor/cmd/...@main`.
  The release installer cannot be used here: the container reaches github.com
  through a proxy that serves only the repositories attached to the session,
  so the release download is refused with HTTP 403. The Go module proxy is
  reachable directly; the image's Go 1.24.7 fetches go1.26.8 by itself.
  `CONDUCTOR_VERSION` overrides `main`.
- **Installs cairn.** It runs
  `cargo install --locked --root /usr/local --git https://github.com/aburan28/cairn cairn`,
  without the `ui` feature, so `cairn mcp` works and `cairn run` does not.
  `CAIRN_REV` pins a commit. `.mcp.json` starts `cairn mcp` for every session,
  on a ledger at `~/.cairn/cairn.jsonl` rather than in the repository; it can
  post only unfunded objectives.
- **Brings up a control plane.** If `CONDUCTOR_ENDPOINT` and `CONDUCTOR_TOKEN`
  are set (see below), the session uses that shared one. Otherwise the hook
  starts a control plane in the container: a Postgres 16 cluster from the
  image on `127.0.0.1:55432`, then `conductor up --project crypto` on
  `http://localhost:8080`. With it the pre-edit hook, `conductor check` and the
  MCP tools work, but it coordinates nothing outside that container.

The first launch is the slow one: 128 s when this was tested, mostly
fetching Go's toolchain and building cairn. The container is cached after the
hook finishes, so later launches skip the installs and only start the
services, which took about a second. A container restored from that cache
keeps a `postmaster.pid` from a server that no longer runs; the hook removes
it when that pid is not a live `postgres`. Everything the hook prints goes to
`~/.local/state/workspace-tools/session-start.log` (mode 0600, since
`conductor up` prints a login token), except one summary line for the
session. Outside the cloud the hook does nothing. On the very first launch,
Conductor's own `conductor hook session-start` entry may run before the binary
exists; it then fails once without blocking anything.

It starts no autonomous agent. `conductor worker` would claim queued tasks and
run Claude on them with `acceptEdits`. `cairn agent run` would execute jobs
other people submit, and this container has none of the sandboxes it uses
(Kata, gVisor or bubblewrap). Start either by hand if you decide you want it.

### Connecting to a shared control plane

A control plane in the container cannot see other sessions. To coordinate
with them, set these as environment variables in the environment's settings,
under the cloud environment menu in the session's title bar, then Edit. The
CLI, the hooks and `conductor-mcp` all read them in place of
`~/.conductor/credentials`.

| Variable | Value |
| --- | --- |
| `CONDUCTOR_ENDPOINT` | the HTTPS URL of a control plane reachable from the internet |
| `CONDUCTOR_TOKEN` | the token of a dedicated account: `conductor member add claude-cloud --role contributor` prints it once, in a `conductor login …` line. It is a credential, so never commit it or paste it into a chat |
| `CONDUCTOR_PROJECT` | `crypto` |

Then add that endpoint's host under **Network access → Allowed domains**,
keeping the package-manager defaults ticked. The endpoint must serve TLS:
`conductord` refuses a reachable address in plaintext (`--tls-cert/--tls-key`,
or `--behind-proxy` behind a proxy that terminates TLS). A loopback daemon on
a laptop, or a Tailscale-only name, cannot be reached from the container.

## Checking it works

```sh
conductor version
conductor doctor
conductor check --summary "what you are about to do" --scope path:docs/CONDUCTOR.md
```

`conductor check` exits 0 when the path is free and 3 when someone else holds
it; the pre-edit hook then blocks the edit and names the holder. Exit 1 means
no control plane answered, or you are not logged in.

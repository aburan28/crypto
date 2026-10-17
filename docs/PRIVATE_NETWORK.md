# Private Conductor and Cairn workspaces

This repository uses [Conductor](https://github.com/aburan28/conductor) for team
coordination and [Cairn](https://github.com/aburan28/cairn) for the research/worker
network. Their native identity, collaboration, verification, and storage layers
remain authoritative. Do not deploy a second account/task system alongside them.

The private path is intended for isolated enterprise installations with approved
local binaries, an internal Conductor/PostgreSQL service, and explicit internal
network endpoints. It makes no network or package installation decision during
workspace startup. A physical air gap is preserved by the organization's network
controls; an application profile cannot enforce a firewall or cross that gap.

## Install and configure before the session

Provision compatible Conductor binaries and a Cairn release with its embedded UI
through the organization's approved offline package process. A compatible Cairn
binary advertises `CAIRN_NETWORK=private` in `cairn run --help`; the launcher rejects
older binaries that could silently ignore the environment setting. Supply the
local verifier bundle and any required model/package artifacts through the same
process. Runtime startup never clones, downloads, builds, or updates them.

Configure these in the environment's secret/configuration mechanism:

```sh
export WORKSPACE_OFFLINE=1
export CAIRN_NETWORK=private
export CONDUCTOR_ENDPOINT=https://conductor.internal:8443
# CONDUCTOR_TOKEN is supplied by the organization's secret mechanism.
export CONDUCTOR_PROJECT=crypto
export CAIRN_PRIVATE_BOOTSTRAP=/etc/cairn/internal-peers.json
export CAIRN_PRIVATE_ROOT=/opt/cairn-verifiers
export CAIRN_PRIVATE_STATE=/var/lib/cairn/workspace
export CAIRN_NODE_P2P=10.88.0.2:9100
export CAIRN_NODE_HTTP=10.88.0.2:9180
```

Only use existing permitted paths and certificates trusted by the internal CA.
Loopback defaults are available for local operation. An internal HTTP Conductor
endpoint is refused unless it is an explicit loopback origin. Tokens are never
written into a tracked file or endpoint URL. The launcher's compatibility check
does not authenticate to Conductor or verify that its service is healthy; use
`conductor doctor`/`conductor status` against the approved internal endpoint.

`CAIRN_PRIVATE_BOOTSTRAP` is an operator-provisioned **local file**, in Cairn's
existing bootstrap format. Do not download it from a seed website. Its identities
form the permitted peer set; configure reciprocal identities on the other nodes.
When omitted, the private node is a standalone leader with no permitted P2P peers.
This is useful for an internal HTTP fleet service, not a discovery mechanism.

## Check and run

```sh
bash scripts/private-workspace.sh check
bash scripts/private-workspace.sh run
```

`check` verifies required local executables, Cairn's embedded reader/profile
support, endpoint configuration, and optional local bootstrap/verifier paths. It
does not create state or start a service. `run` creates a mode-0700 state directory
and execs the native Cairn process in the foreground. Use your approved service
manager for supervision, logs, process ownership, and restart policy.

The launcher forces `CAIRN_NETWORK=private`, `CAIRN_SEEDS=off`, bootstrap-only peer
policy, and disables beacons and port mapping. The native Cairn profile additionally
suppresses public seed refresh and ambient tailnet discovery. It retains native
member/node authentication and rejects incompatible public/auth settings rather
than silently removing authentication. Set `CAIRN_NETWORK=private` for other
independently launched Cairn processes too; a shell child cannot set its parent's
environment. The launcher starts no autonomous worker or spending agent.
It also removes the Conductor bearer credential from Cairn's environment,
including its local feature probes; the caller's own session keeps that credential.

With `WORKSPACE_OFFLINE=1`, `.claude/hooks/session-start.sh`,
`scripts/setup-cloud-environment.sh`, and `scripts/setup-workstation.sh` run the
local check and return before their online installation/startup paths. Missing or
incompatible tools fail with an actionable message instead of falling back to
public installation. Offline mode also avoids automatic PostgreSQL creation;
provision the approved internal service separately. Unset or `0` preserves the
existing online setup behavior. Invalid values fail before any installation.

The existing public `scripts/cairn-node.sh` and its mirrored copies in sibling
repositories are unchanged. Do not invoke that public launcher for an isolated
installation. Other explicitly invoked online installers and tools remain the
operator's responsibility; `WORKSPACE_OFFLINE` is not an OS-wide network switch.

## WireGuard and service discovery

WireGuard supplies encrypted private reachability. Conductor/Cairn retain their
application authentication and authorization. Generate an operator-reviewed plan
using the native Conductor renderer:

```sh
conductor wireguard render --manifest /etc/conductor/mesh.yaml --out ./wireguard-plan
```

The strict version-1 manifest describes a private overlay CIDR, one gateway's
public key/address/endpoint, and named node public keys/addresses. It contains no
private keys. The renderer creates owner-only gateway/node configuration templates
in a new directory, refuses overwrite, uses exact /32 gateway peer routes, and
routes workers only to the private overlay. Replace private-key placeholders
locally on each machine and apply configurations through approved host tools.
See Conductor's `docs/WIREGUARD.md` for the complete manifest and limitations.

The gateway requires an explicit host firewall and forwarding policy before use.
Templates do not provide workspace-level packet isolation or internet exit/subnet
routing. A WireGuard gateway decrypts and re-encrypts forwarded packets; use
application TLS when payloads must remain private from its operator. Revoking an
application identity does not remove a tunnel peer; network administrators must
reconcile gateway membership too.

Use internal DNS or static addresses for controller and peer names. The native
mesh DNS can be configured with an explicit overlay address and hostname, avoiding
Tailscale discovery. Its SRV target must resolve through the caller's internal
system resolver; configuring a custom mesh SRV resolver does not automatically
replace HTTPS resolution. See Conductor's `docs/AIR_GAPPED.md`. Existing mesh
discovery establishes reachability; it does not federate organizations or tasks
between independent Conductor databases.

## Verification and remaining enterprise work

### Repeatable startup acceptance

The fixture contract suite runs on every matching PR in
`.github/workflows/private-workspace.yml`. It includes 32 rejection cases plus
successful launches and exact argument checks. Run it with:

```sh
bash scripts/test-private-workspace.sh
```

For a Linux release bundle already installed on `PATH`, run the separate real
process acceptance test:

```sh
PATH=/opt/approved-cairn/bin:/opt/approved-conductor/bin:$PATH \
  bash scripts/test-private-workspace-native.sh
```

This test requires Python 3 and readable Linux `/proc`. It installs nothing,
uses temporary state and loopback ports, records the tested Cairn binary digest,
starts and stops only its own process group, and removes its state afterward.
It uses a dummy loopback Conductor origin and sentinel credential: it does **not**
authenticate to a control plane. Native acceptance is opt-in; the GitHub workflow
runs fixture contracts and shell syntax, not a cross-repository release build.

| Behavior | Executable check | Boundary |
| --- | --- | --- |
| All three offline entrypoints avoid installers and service/auth startup | Fixture tripwires and no state/process artifact | Commands intercepted by the fixture, not an OS firewall |
| Missing/old/headless/failed capability probes refuse startup | Fixture rejection cases, including nonzero probes that print valid markers | Approved real release is checked separately |
| Origins, ports, bootstrap and verifier paths fail closed | Fixture invalid/boundary inputs; no launch/state before rejection | Local syntax/path checks, not service health/authentication |
| Public environment cannot weaken launch | Fixture child assertions; real `/network` and empty seeds | Does not prove physical network isolation |
| Least authority, exact arguments, permissions and exit status | Fixture relay/spend limits, default/custom argv, reused state mode 0700, native failure exit 37 | Host service-manager policy remains external |
| Bearer credential isolation | Both fixture probes/child and real `/proc` assertion | Caller credential remains available to Conductor |
| Packaged reader, private health and identity survive restart | Two real launches, HTML and local JavaScript fetches, identity equality | Same Linux host and supplied binary; browser interaction covered in Cairn |
| Real startup failures reach the supervisor | Occupied HTTP port causes nonzero exit with bind diagnostic; identity preserved | Tests failure, not an automatic recovery service |

These are behavior coverage obligations, not a whole-project line-coverage
percentage. Conductor's native WireGuard/team tests and Cairn's peer/browser
tests live in their respective repositories; passing shell contracts alone does
not establish either of those layers.

Local Linux validation on 2026-10-17 passed shell syntax, the 32 rejection cases
and successful fixture launches, and the complete native acceptance script on the
fresh Cairn reader/recovery candidate (binary SHA-256
`0c24a34d932fd58a9760270a4129c31eadbfd9b32a3929a3ca0ba021ed3d2888`).
The unchanged Rust library at base `b952011cfd9a871c2861d8d11d2c5c3fae65844a`
also passed `cargo test --release --lib -- --test-threads=1`: 4,268 passed,
zero failed, 117 existing ignored, 926.77 seconds. These are local results;
GitHub checks and customer acceptance remain separate gates.

Run `bash scripts/test-private-workspace.sh`. It exercises all three entrypoints
with install/network commands replaced by failing tripwires, rejects legacy
profile support and invalid configuration, and checks the private environment and
explicit bootstrap arguments passed to the native node. The new CI workflow runs
the same checks. A real embedded Cairn binary also passed the launcher's local
health and empty-seed checks. Native Conductor-rendered configurations established
two authenticated WireGuard gateway handshakes and delivered 12/12 test packets
through an isolated Linux userspace hub/worker network; explicitly removing a
gateway peer then made it unreachable. The first forwarding attempt exposed ICMP
redirect packet loss, so the deployment instructions include redirect controls.
See Conductor's
[routing validation and preserved evidence](https://github.com/aburan28/conductor/blob/5c4f63cc08570ecb5c1a8dcb428ee8caf30f64de/docs/WIREGUARD_VALIDATION.md).
These checks do not establish kernel WireGuard compatibility, customer firewall
isolation, automatic revocation propagation, or production availability.

| Requested foundation | Status and boundary |
| --- | --- |
| Shared team work, role isolation, identity | Reuse existing Conductor/Cairn implementations; no replacement identity database introduced |
| Private workspace startup | Implemented; local shell contracts verify no installer/public-startup fallback |
| Private peer discovery | Native Cairn private profile with explicit bootstrap; no automatic public rendezvous |
| WireGuard routing nodes | Native static configuration renderer; host activation and firewall remain operator-owned |
| Internal service discovery | Explicit native mesh DNS and existing Cairn bootstrap/registry; requires internal name resolution |
| Physically isolated deployment | Requires customer OS/network acceptance with outbound denial and supplied CA/time/model/package services |
| Automatic mesh, gateways/exits, federation | Not implemented by this change; do not infer them from private reachability |

Production acceptance must repeat the local routing, authorization, restart and
revocation checks on the intended customer topology. It must also exercise kernel
gateways, denied cross-workspace packet routes, key rotation, signed offline
updates, customer identity/secret-provider configuration, and workload scale.
Private mode is not a claim of production readiness or additional cryptographic
security. Application authorization tests do not establish firewall isolation
between workspaces or automatic reconciliation of WireGuard membership.

Coordination during development: `conductor check` was attempted before edits;
the CLI was initially unavailable. Work proceeded on an isolated feature branch.
After a native CLI became available, its check was retried; the resulting access
status is recorded in the PR. No chat transcript, credential, or prompt is task
metadata.

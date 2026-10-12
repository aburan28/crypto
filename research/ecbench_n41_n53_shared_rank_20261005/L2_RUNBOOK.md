# L2 host runbook for the frozen n37 and n41/n53 shared-rank specs

This runbook carries two frozen `ecbench.spec/v1` files to an auditable
physically isolated host so that their **wall-clock** figures can be read at
the level the specs require (`isolation_required: L2`). Operation counts do
not need it; they are already settled by the macOS L0 sessions and replayed
by CI. What the L2 run adds is the primary one-target online wall comparison
(AGENTS.md "IC measurements"; `docs/ecbench/README.md` §3 and §6), which
neither a GitHub-hosted VM (L1, `research/ecbench_n37_native_online_wall_20261004`)
nor this Mac (L0, `sessions/mac-l0`) can supply.

## What runs

| Spec | Unchanged file | Curves | Executions |
|:--|:--|:--|--:|
| n37 native online wall screen | `research/ecbench_n37_native_online_wall_20261004/SPEC.json` | `icv1-f2m37-tm534059-32aad96b` | 384 |
| n41/n53 shared-rank counted panel | `research/ecbench_n41_n53_shared_rank_20261005/SPEC.json` | `icv1-f2m41-tm2308219-7f48b14a`, `icv1-f2m53-tm56619371-dac20a85` | 768 |

Both specs are run **byte for byte as committed**; the host run changes no
seed, cap, arm or workload. The n37 spec is the one whose decision was
`carry_both_host_noise_exceeds_gate`; the n41/n53 spec is this round's.

## Host

- Linux bare metal, Ubuntu 24.04, run as root (core reservation and thread
  eviction need it; a non-root run stops at L0).
- The intended host is an AWS bare-metal instance in `us-west-2` launched with
  the existing key pair **`meow34`** (AGENTS.md §9: `--key-name meow34`, SSH
  with `meow34.pem` at mode 0600, never commit or print the key).
- **Launching that host needs the repository owner's authorization and was
  not performed in this round.** This round was explicitly denied any cloud
  launch or paid resource; everything measured here ran on the Mac at L0.
- A VM is not a substitute: a hypervisor cannot certify zero steal, and the
  n37 hosted run already showed what an L1 VM yields.

## Procedure

`run-on-host.sh` does the following, in order, and stops on the first error:

1. Installs the build toolchain and `rustup` (minimal profile) if `cargo` is
   absent.
2. Clones `aburan28/crypto` and checks out the exact commit named in
   `COMMIT` (the merge commit of this round's PR, or any later commit on
   `main` that leaves both specs and the `ecbench` method ids unchanged).
3. Builds `ecbench` in release mode once and records the binary SHA-256 and
   `rustc --version`; every record carries that binary hash.
4. Quiets the host: stops and disables `unattended-upgrades`, the `apt-daily`
   timers and services, `snapd`, `fwupd-refresh`, `motd-news`, `man-db` and
   `e2scrub_all` timers; sets the `performance` governor on every CPU; turns
   turbo/boost off where the sysfs knob exists; disables swap and NUMA
   balancing. Each step is best effort and its output is kept in
   `quiet-<stamp>.log`; `ecbench`'s own preflight and per-run grading decide
   the level each run actually earns.
5. Records `ecbench host`, `uname`, `lscpu`, `/proc/cmdline` and the governor
   beside the sessions.
6. Runs `ecbench plan --json` for both specs, then `ecbench run` for each with
   `--cpus auto --wait` and **no `--allow-busy`**: a busy preflight refuses to
   start, which is the gate working.
7. Runs `ecbench verify --replay-all --exit-code` on both sessions and writes
   the receipts (`n37-AUDIT.json`, `n41-n53-AUDIT.json`), then `ecbench table`.
8. Hashes every file into `SHA256SUMS` and tars the session directory.

```sh
# on the host, as root
COMMIT=<merge commit of this round> ./run-on-host.sh /root/ecbench-l2
```

Expected duration on an idle bare-metal core: the n37 spec is a few minutes;
the n41/n53 spec is dominated by its IC arms. On the busy Mac at L0 the
n41/n53 panel needed the wall recorded in `RESULT.md`; budget at least that
on the host, and do not shrink `timeout_seconds` to fit.

## What to bring back, and what it is allowed to say

- Copy `ecbench-l2-sessions-<stamp>.tar.gz` off the host, record its SHA-256,
  and commit the two session directories under
  `research/ecbench_n37_native_online_wall_20261004/sessions/` and
  `research/ecbench_n41_n53_shared_rank_20261005/sessions/` (new
  subdirectories; never overwrite `hosted_ubuntu_01` or `mac-l0`).
- Replay both sessions on a different host class (`ecbench verify
  --replay-all` on the Mac or in CI) and commit those receipts too; the
  receipt's SHA-256 is the replay certificate a claim cites.
- Re-run the frozen analyzers. For n37:
  `examples/ecbench_n37_native_online_wall_analyze.rs` with the new session
  and its independent receipt, under the gate its protocol froze (all
  measured pairs at L2, K16 A/A maximum paired deviation below 5%). For
  n41/n53: `examples/ecbench_n41_n53_shared_rank_analyze.rs`, which reports
  the counted quotients (host independent) and the online-wall ratios with
  the levels the host actually earned.
- A wall figure is admitted only if every pair it rests on earned L2 and the
  A/A gate passed; otherwise it stays a descriptive diagnostic, exactly as
  the n37 hosted run did. The counted cold IC/rho lower-bound quotients are
  the same on every host and are not a wall result.
- Any promotion goes through the scoreboard rules of AGENTS.md §7 and §7a in
  the PR that lands the sessions.

Status of this runbook: **pending; the host has not been launched.**

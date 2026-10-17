# Backend leaderboard reproduction contract

An independent agent must identify **which result** it reproduced. Use the
three gates below separately. A regenerated page proves nothing about a new
solver execution. Matching operation counts proves nothing about milliseconds
on different hardware. A missing or failed gate stays open.

## Acceptance criteria

| Gate | Required agreement | Pass evidence | Failure action |
|---|---|---|---|
| R0: report reproduction | Every pinned input and all five generated files match SHA-256 exactly; summary and frontier memberships agree | Exit 0 from `reproduce` using a lock hash obtained from a trusted reviewed commit | Stop; keep diagnostics. Do not reseal to make a reproduction pass |
| R1: arithmetic and deterministic ledger reproduction | Same exact spec, plan, target, curve, field representation, subgroup, methods, parameters, seeds, budgets and attempt count. All measurements verify. Recovered scalars, charged totals, exclusive phase counts, native counters, factor-base identity, field operations when present, and both recorded solver-statistics blocks agree exactly, excluding explicitly listed clocks | Fresh native `--replay-all` audits of original and new sessions, actual executable hash, plus exit 0 from `check-rerun` | Retain the discrepancy, timeout, OOM or failed check. No averaging, dropping attempts, missing-to-zero conversion or relaxed epsilon |
| R2: timing reproduction | Same declared environment class and timing window, complete L2-or-higher observations, preregistered A/A controls and paired equivalence test below | New raw sessions, host/build/placement evidence and statistical receipt meeting every timing criterion | Report descriptive times or failed/inconclusive replication; no timing badge |

R0 is automated in CI. `check-rerun` automates R1 **validation of supplied
records and audit assertions**. It does not run the solver itself, authenticate
someone else's receipt, prove physical host independence, or admit R2. The
external runner's job log, commands and retained executable are the execution
witnesses. Hashes detect changes relative to a trusted commit; someone who
rewrites the code, lock and receipts can forge a self-consistent story. Obtain
the commit/hash through the reviewed repository, not the submitting agent's
working tree alone.

## R0: exact report reproduction

Check out the full 40-character reviewed PR/merge commit, not moving `main`.
The lock is `docs/backend-leaderboard/reproduction.lock.json`. It pins the
manifest, all five files of each imported session, each historical receipt,
every producer source/template/build input, Cargo dependencies, Rust toolchain,
the expected result counts and the hashes of:

- `leaderboard.json`
- `LEADERBOARD.md`
- `index.html`
- `REPLAY_PLAN.json`
- `REPRODUCING.md` (the published copy of this contract)

Materialize paths from the manifest in sparse checkouts. A missing blob is an
error, not permission to skip a source. From the repository root, substitute
`TRUSTED_LOCK_SHA256` with the SHA-256 of the lock in that reviewed commit:

```sh
cargo +1.99.1 fetch --locked --manifest-path tools/backend-leaderboard/Cargo.toml
cargo +1.99.1 test --locked --offline --manifest-path tools/backend-leaderboard/Cargo.toml
cargo +1.99.1 run --locked --offline --manifest-path tools/backend-leaderboard/Cargo.toml -- \
  reproduce --expected-lock-sha256 TRUSTED_LOCK_SHA256 > report-reproduction.json
cargo +1.99.1 run --locked --offline --manifest-path tools/backend-leaderboard/Cargo.toml -- --check
```

`reproduce` recomputes every view in memory and checks pinned output hashes; it
does not trust the checked-in generated views. `--check` additionally checks
those views on disk and never repairs them. Both fail closed for altered inputs,
an incomplete lock, a stale executable, a changed compiler or changed results.
The receipt says `solver_executed: false` and `fresh_arithmetic_replay: false`.
Exact report bytes are required across hosts; compiler executables and PDF
metadata are not claimed to be reproducible builds.

For a deliberate new baseline, rebuild, run `--seal`, review the input and
expected-result changes, and commit all five generated files and the lock together.
`--seal` is a maintainer action, never a replication repair step. CI uses the
already committed lock and cannot silently update it. Cargo's lockfile checksums
pin dependencies; vendor the locked crates in a separately hash-bound bundle
when the reproduction machine must begin fully offline.

## R1: rerun the algorithms

`REPLAY_PLAN.json` contains one command card per exact imported session, its
source/binary identity, all five input hashes, the exact required replay
record IDs and known prerequisites. Command cards are **planned instructions**,
not evidence that a fresh execution occurred. The acceptance unit of
`check-rerun` is an entire session; a strict pass cannot hide another arm's
failed or nondeterministic attempt. Preserve blocked cells as blocked. A new
arm-only experiment needs its own preregistered spec and cannot claim to have
reproduced the full old session.

Before executing, freeze and retain:

1. The reviewed report commit and lock SHA-256; exact session/spec/plan hashes.
2. The original executable, **or** the complete clean source revision and
   dependency lock/snapshot needed to rebuild it. Record `rustc -Vv`,
   `cargo -V`, target triple, dependency hashes, features, compiler/linker flags
   and the build command with its output/exit status. The reporting tool's
   Rust 1.99.1 pin does not describe how historical solver binaries were built.
3. The actual `ecbench` executable and SHA-256, all external solver versions
   and executable/library hashes, and `ecbench host --json` output. For a
   rebuilt binary, both host build capsules and session declarations must
   identify the same clean **build-time** source commit, with
   `git_commit_source` equal to `CRYPTO_BUILD_GIT_COMMIT` or `GITHUB_SHA`.
   `runtime_worktree`, null or missing provenance is insufficient: a checkout
   beside a prebuilt executable does not identify its source. A dirty or
   unknown historical tree requires its exact
   original binary or a complete hash-bound snapshot; a commit alone is
   insufficient. The checker currently supports exact binary or matching
   clean build-time commit; resolving a dirty snapshot requires a separately reviewed
   new protocol, not silently changing this gate. Set build provenance only
   after independently checking the actual clean source and pinned build inputs;
   declaring a revision in an environment variable is not a substitute for that
   check. Retain the build/job log that witnesses it.
4. Seeds, SAT variable/branching/restart/XOR/conflict/learned-clause settings,
   F4/F5/F6 and root-finding configurations, factor base, linear algebra,
   operation and memory budgets, thread counts, CPU affinity, NUMA, warmups,
   ordering and timeout policy. The exact spec and plan carry these inputs;
   no command-line override may quietly change them. Pin undocumented external
   randomness before calling its counters deterministic.
5. A new evidence directory on a separate runner. Preserve argv, cwd, start/end
   times, exit codes, stdout/stderr, host capsule and job/artifact location for
   every command, including failed commands. Never overwrite the original.

For example, select the K13 session. Run from the repository root, with `ecbench`
already built from the frozen source and `rerun-evidence/` a **new** directory:

```sh
mkdir rerun-evidence
ecbench host --json > rerun-evidence/host.json
ecbench run --spec research/f4_f5_f6_cross_family_20261011/sessions/K13/spec.json \
  --out rerun-evidence/session
ecbench verify --dir research/f4_f5_f6_cross_family_20261011/sessions/K13 \
  --replay-all --exit-code --out rerun-evidence/original.audit.json
ecbench verify --dir rerun-evidence/session --replay-all --exit-code \
  --out rerun-evidence/rerun.audit.json
cargo +1.99.1 run --locked --offline --manifest-path tools/backend-leaderboard/Cargo.toml -- \
  check-rerun --source research/f4_f5_f6_cross_family_20261011/sessions/K13 \
  --rerun rerun-evidence/session \
  --original-audit rerun-evidence/original.audit.json \
  --rerun-audit rerun-evidence/rerun.audit.json --binary /absolute/path/to/ecbench \
  --expected-lock-sha256 TRUSTED_LOCK_SHA256 > rerun-evidence/ledger-reproduction.json
```

Use the same absolute executable for all three `ecbench` invocations. Every
native audit must exit 0, bind all five files, have no failed checks, identify
that executable, and reproduce **exactly the entire expected set**, with no
duplicates. Both audits must be generated after the fresh session finishes;
the new session must start after the old session finished. Legacy audit formats
without `records_bound` must explicitly pass `records_hash`, in addition to
the direct binding checks. `--replay 3` or a successful subset does not qualify. The checker
also rechecks file/seal/plan integrity and compares fresh raw records; the native
audit's old comparison contract does not cover the complete solver block.

Ledger equality excludes session/run/record IDs, timestamps, placement, RSS
and clocks. Among each phase's direct fields and each solver block's direct
fields, only keys named `wall_ns`,
`wall_seconds`, `wall_ms`, `cpu_ns`, `solve_wall_ns`, `elapsed_ns`,
`elapsed_ms` and `elapsed_seconds` are excluded. Nested fields and counters
are retained even if they have one of those names. No counter is excluded.
Both top-level `solver` and `detail.decomposition.solver` are compared, so
conflicts, learned clauses, decisions, propagations, work units, solver budgets
and reported outcomes cannot silently change. A newly missing observation fails.
Other method-specific detail is retained in raw files but is **not** asserted
to be covered by ledger equality; adding a new claimed metric requires extending
and versioning this contract and its tests first. Warmup ledgers also compare,
although native `--replay-all` selects measured deterministic verified runs.

Counts use exact integers and JSON float representations, with no tolerance.
Historical arithmetic audits may describe their explicit legacy rounding rule;
the fresh record comparison still requires equal recorded costs. A new counting
rule or instrumentation schema needs a new method/protocol and evidence.

An R1 pass is a deterministic-ledger reproduction. To call it **independent**,
also retain a separately controlled runner/job identity and durable artifacts,
and compare original/auditor hardware profiles. Different environment IDs alone
can reflect a compiler or OS change on the same machine; they do not prove a
second host. No independent-host badge is inferred by this tool.

## R2: reproduce time or beat a timing frontier

Preregister the following policy before collecting new timing data. It is a
stricter replication design than the report's minimum five-pair historical
timing-admission check; it does not retroactively qualify old observations.

- Use the same environment class, architecture, CPU/microcode/frequency/SMT,
  RAM/NUMA, OS, compiler policy, thread count and resource envelope. Record
  actual earned isolation levels and blockers, not only requested placement.
- For every exact target/configuration, collect at least **20 interleaved
  matched rounds** and at least 20 same-target same-method A/A pairs at L2 or
  higher. Fix the target/algorithm seeds, warmup policy and randomized order in
  the preregistration. Keep every attempt. A timeout, OOM, failed answer or
  incomplete phase window blocks that comparison.
- Declare **5% relative equivalence** in advance. Compute each paired
  `rerun_time / original_time` ratio for the declared window, then its median
  and paired 95% percentile-bootstrap interval with 10,000 resamples and fixed
  seed 20261017. The entire interval must lie within **[0.95, 1.05]** for a
  timing reproduction. Merely failing to reject equality is insufficient.
  The A/A interval on each campaign must also fit inside that band; otherwise
  improve control and start a newly documented campaign. Never widen the band
  after seeing the result. A different tolerance needs a new preregistration.
- For an improvement, retain the repository's same-target paired interval
  excluding parity **and** the A/A noise interval. Require the preregistered
  full-operation improvement on untouched holdouts, not just fewer conflicts,
  clauses or matrix rows. Publish every target's result and multiplicity policy
  before making a many-target/general claim.
- Keep cold setup-to-recovery, fork-to-reap, online and independent verification
  boundaries separate. Include reusable setup/construction when claiming the
  full cold operation. Report excluded independent scalar-check cost. Never
  substitute online-only time or a batch average for a single-target total.
- Compare the fastest available external tool on identical input and host
  before a fastest-tool claim. Record missing Magma/F4/F5, msolve, FGb or
  FLINT/NTL/root references explicitly. Missing references block that claim.

This tool does not issue an R2 verdict. Save the statistical implementation,
paired data, intervals and manual gate review as a separate receipt. Current
imported sessions have no admitted timing winner and do not satisfy this R2
protocol. Different-host milliseconds are descriptive and cannot satisfy it.

## Current completion and stop conditions

R0 is implemented and tested against the committed panel. R1 comparison and
adversarial checks are implemented; **no fresh algorithm campaign or independent
host rerun is asserted by this change**. R2 remains unmeasured. The current
root integration cannot build `ecbench`; see the
[tracked compiler blocker](../../problems/ecbench-main-build-integration.md).
Missing historical source/binary dependencies must be recovered before their
command cards can be executed. None of the 17 imported host capsules identifies
a clean build-time source under this strict gate: the observed provenance is
`runtime_worktree` or absent. Six ladder sessions also name no source commit,
and the degree-19 root session declares a dirty tree. Their exact original
executables are required by the current R1 checker. Do not
use the report's reproducible bytes as evidence those binaries can be rebuilt.

Keep a requirement-to-evidence checklist for R0, R1, independent-runner proof,
R2, external references, holdouts and every requested field size. Degrees
3/4/5/51/83 have no imported measurements. No run at a smaller degree discharges
those cells or the exact E₀ degree-83 gate. Stop on a build error, a mismatched
hash/counter, failed correctness, exceeded preregistered budget, or inadequate
timing control; preserve outputs and report the precise open gate.

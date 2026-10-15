# cairn certificate objectives pin a Python checker; crypto forbids Python verifiers

- **Status:** `open`
- **Kind:** implementation gap (in cairn, a separate repository); policy conflict
  between two repositories' rules. Not a mathematical obstruction.
- **Affected work:** objective mode for any backlog card with a native
  checker: the upgrade path of
  [`backlog/XP-20261015-p192-isogeny-aux-search-t48`](../backlog/XP-20261015-p192-isogeny-aux-search-t48/README.md),
  which since 2026-10-15 runs as cairn lab tasks with the native checker
  deciding `found` and is no longer blocked by this record.
- **Next action owner:** unassigned (the change, if any, is in cairn; this
  repository cannot push there).
- **Last updated:** 2026-10-15

## Problem

We cannot post the P-192 auxiliary-curve search as a cairn `certificate`
objective with a native (Rust) checker, because cairn runs every pinned
certificate checker as a Python module, and AGENTS.md ("Implementation
language: no Python") forbids Python verifiers in this repository.

The backlog's design (backlog/README.md) puts checkable searches on cairn
objectives: a contributor runs `cairn work` with any solver, and a pinned
verifier, not the submitter, decides whether a claimed curve is valid. The
check for this search is cheap and fully specified (artifact `{a, b, Q, t}`;
`p` is the P-192 prime, `4a³ + 27b² ≠ 0`, `Q` on the curve, `[2^t]Q = O`,
`[2^{t−1}]Q ≠ O`), so the only gap is where the checker may run.

## Evidence and uncertainty

- **Observed (cairn source, read-only clone at `aburan28/cairn`, the version
  installed as `cairn` 8.0.0):**
  - `src/scaffold.rs:288–296`: a `certificate` verifier block is
    `{kind, checker: <path>, checker_sha256, entrypoint}`; the pin is the
    SHA-256 of a script body.
  - `src/verifiers/mod.rs:57–60` (module docs): "The pinned checkers in
    `examples/` are Python. Rust cannot `exec` them; it can spawn an
    interpreter that can."
  - `src/verifiers/mod.rs:916–1000`: the checker runs as `python3 -c
    HARNESS <pinned file> <entrypoint>`; `HARNESS` loads the pinned file with
    `importlib`'s `SourceFileLoader`, calls the entrypoint on the artifact
    read from stdin, and reads `{ok, detail}` back. There is no branch for an
    executable.
  - `src/verifiers/mod.rs:1118–1150`: the registry has `python_binary` and
    `lean_binary` and no native-checker setting.
  - Every certificate objective under `examples/` pins a `.py` checker (for
    example `examples/certicom-ecdlp/objective-ecc2k130-orbit-batch.json`).
- **Derived:** the two other verifier kinds that run commands do not close the
  gap as they stand.
  - `replay` (`src/verifiers/mod.rs:2261–2420`) runs a command resolved on the
    operator's `PATH` (not pinned by hash) with no stdin, and compares fields
    of its stdout with `artifact.results`. The command never sees the
    submitted witness, so it cannot check a submitted curve.
  - `workspace` (`src/verifiers/workspace.rs`) applies the artifact's files
    over a pinned base tree (`base_sha256`) and runs `score_command` with the
    tree as its working directory, jailed and without network; the score file
    must equal `artifact.results.score`. A base tree holding a zero-dependency
    Rust checker crate, with `score_command = ["cargo", "run", "--offline",
    "--release", "-q", "--", …]`, would pin the checker's *source* and run it
    natively. The program (`cargo`) is resolved on the operator's `PATH`
    (`workspace.rs:388`), as `lean` is for Lean objectives.
- **Not established:** whether `cargo` and `rustc` start inside cairn's jail
  (gVisor or bubblewrap, read-only host binds of the resolved binary and the
  bundle root only). A rustup-managed toolchain reads `~/.rustup`, which the
  jail does not bind, so the likely outcome is `unavailable`. This host has
  neither `runsc` nor `bwrap` (`cairn lab sandbox`: "no sandbox works here"),
  so it could not be tried here.
- **What this does not establish:** nothing about the search itself, its cost
  or its outcome.

## What would unblock it

Either of:

- [ ] **cairn accepts a native pinned checker** for `certificate` objectives:
  the objective pins an executable (or a static binary blob) by SHA-256, the
  verifier runs it jailed with the artifact on stdin and reads an exit code or
  `{ok, detail}` from stdout. Acceptance: a certificate objective whose
  checker is a static Rust binary settles an accepted and a rejected claim on
  two nodes, and a tampered binary yields `invalid_spec`. This is a change to
  cairn and the user's decision.
- [ ] **The `workspace` route works in the jail:** a base tree with the Rust
  checker crate (no dependencies, so `--offline` builds) and `cargo` as
  `score_command` reaches `accept` for a valid curve and `reject` for an
  invalid one on a host with gVisor or bubblewrap. Acceptance: those two
  verdicts plus `unavailable` (never `reject`) on a host without a Rust
  toolchain. Limitation that would remain: the toolchain is the operator's
  and unpinned, as Lean's is.

Meeting either lets the card's `cairn.checker.status` become `available`. The
search stays blocked on its own prerequisites (see the card).

## Next discriminating step

- **Action:** on a host with `bwrap` or `runsc`, scaffold a `workspace`
  objective with a ten-line Rust crate whose `main` writes `{"score": 1}` and
  run `cairn propose <objective> --artifact … --dry-run`.
- **Outcomes:** `accept` → the workspace route is viable; write the real
  checker crate (no dependencies) and pin it. `unavailable` naming the
  toolchain → only a cairn change unblocks it; raise it with the cairn owner.
- **Recheck trigger:** a cairn release that mentions native checkers, or the
  P-192 card's other blockers clearing.

A planned step is not a completed experiment.

## Resolution and history

- 2026-10-15: recorded by the backlog lane while integrating the experiment
  backlog with cairn. Unresolved.
- 2026-10-15 (later): the native checker exists
  (`src/bin/backlog_check_aux2t.rs`, tested on small-prime analogues), and the
  P-192 card switched to `mode: lab-task` with the first-success stop rule of
  `backlog/README.md`, so the search no longer waits on this record. What
  stays blocked is posting it as a cairn objective, where a stranger's claim
  is graded by a verifier the objective pins. Whether to extend cairn so a
  certificate verifier can pin a native executable is listed among the
  user's decisions in `backlog/README.md`.

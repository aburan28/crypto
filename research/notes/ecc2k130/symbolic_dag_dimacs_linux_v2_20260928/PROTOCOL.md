# Held Linux v2 single-edge DAG-to-DIMACS attempt

**Status: HELD. No measured attempt-2 child has been dispatched.** This
successor starts from merged [#804](https://github.com/aburan28/crypto/pull/804)
and [#831](https://github.com/aburan28/crypto/pull/831). #804's immutable first
archive records a macOS `RLIMIT_AS` `LAUNCH_ERROR` before `produce.py` executed;
it is not an algorithmic negative. #831's one Ubuntu 24.04 preparation run
established hard address-space caps and reproducibly built CaDiCaL 3.0.1, but
did not launch toy, n13 or n131. This PR freezes a separate attempt-2 wrapper,
supervisor and replay. `FROZEN.json` deliberately has `status=HELD`,
`release_main_head=null` and `release_pr_number=null`; its measured entrypoint
fails closed until a separately reviewed release commit.

The narrow hypothesis is that the *unchanged* #804 single-edge semantics can
complete on the prepared Linux host: exhaustive toy CNF equivalence, then the
fixed 20-case n13 local group-law SAT/UNSAT panel, then a capped n131
**single-edge export-only** smoke. A pass establishes only those stage gates.
It does not solve a PDP, train a relation matrix, recover a logarithm, measure
an index-calculus total cost `S`, or beat matched rho. A cap, launch error,
missing proof, partial model or incomplete archive is censored or failed,
never a no-go proof for index calculus.

## Frozen identities and only allowed host correction

`FROZEN.json` pins #804's final PR head/merge commit, original first-run
`FROZEN.json`, raw failure receipt and manifest, all original source hashes,
`INPUT.json`, and #831's final head/merge, `PREPARE.json`, preparation freeze,
18-file preparation manifest, hard-cap/build receipts, official CaDiCaL source
commit/tree/manifest and retained executable. The Linux executable is the
committed mode-100755 Git blob
`3f50d3a4218609874b01e02bdf61ad699e6b3f71`, 1,732,064 bytes, with
SHA-256 `dc5dc5754a25e03c3459d2e3220b6b5cc651b7b1628cdbacf76ee5ae12c57a51`.
`ci_replay.py` checks those bytes and invokes #831's independent preparation
verifier. Changing any semantic input or source is a new experiment.

#804 `produce.py` has a latent, masked entrypoint error: its `main()` calls
`release_gate(frozen)` although that function requires an expected-head
argument. The v2 `v2_child.py` imports #804's `toy`, `panel` and `n131`
functions and calls them directly; it never calls the defective `main()`.
Its sole phase-configuration override changes the original Homebrew
`cadical_path` to #831's pinned Linux binary. The original `INPUT.json`,
`bounded.py`, `export.py`, `produce.py`, `verify.py`, DRAT-trim source and all
caps remain byte-identical. `run.py` supplies a new release gate, raw archive
and failure-preserving supervisor. The exact new wrapper, supervisor,
replay, harmless held cap control, protocol and both workflows are SHA-pinned
by `FROZEN.json`.

## Release and one-shot gate

The hash-only workflow checks out the PR's exact head, checks all frozen bytes,
compiles the v2 scripts, runs nine no-network one-shot/refusal and
relation-mutation controls, replays any committed archive, and executes *both*
`v2_child.py --gate-only` and `run.py --gate-only --expected-head <event head>`.
Neither gate-only command executes a measured phase. The outer supervisor
runs full Git ancestry, GitHub API and preparation/archive checks before
dispatch, outside any child address-space cap. Its sealed `DISPATCH.json`
records the reviewed head, event, one-shot run, probe and frozen bytes.
Each capped child consumes that exact dispatch digest, verifies its local
source/input/solver bytes, current checkout head and run identity, then starts
one phase. It cannot start while `status=HELD`. This split is necessary:
on the Ubuntu target image, `git merge-base` exits 128 under the toy
512-MiB `RLIMIT_AS` because Git cannot map the checkout packfile. Repeating
the full gate inside the toy child would falsely refuse a valid ancestor.

The hash-only Ubuntu workflow exercises the actual Git-free child byte gate
under that hard toy cap. It also probes Git, GitHub PR/Actions API and the
preparation gate under the same cap, recording the pack-map refusal as a
specific host limitation rather than an ancestry result. Its harmless receipt
archives command outcomes and zero measured children. A new or different
failure refuses the held control. The separately opt-in measure workflow is triggered only by the unique PR label
`ecc2k130-dag-linux-v2-measure-once`. It cannot release this held freeze.

Before that label can be applied, a new commit must fill `status=RELEASED`,
`release_main_head` with the exact main SHA, and `release_pr_number` with this
PR's number; that new exact PR head, complete diff and hash-only CI must be
independently reviewed. The PR must be ready for review. At dispatch the
runner rechecks the exact checkout/event/live PR heads, #804/#831 merge
identities, Linux x86-64/Python 3.12, the executable, monotonic main ancestry,
and all guarded paths. It rejects any guarded commit between release and live
main, including a changed-then-reverted path. It requires the current Actions
`run_attempt=1` and queries *all* runs of the dedicated label-only workflow by the
current run’s numeric workflow ID:
the current run must be its sole run for this PR branch. The branch check is
necessary because the successful #831 labeled Actions run reports an empty
`pull_requests` array in GitHub’s run API. Even a failed first
label event blocks another label event on the same PR. Other labels on this
workflow are conservatively blocking too; a new reviewed PR/release is needed
if the one-shot event is consumed without a usable archive. A refusal is
recorded before any child and uploaded by the workflow. If checkout, Python
setup or frozen-byte preflight fails before the supervisor can start, an
`always()` workflow step writes a separate pre-dispatch refusal JSON with the
exact run/head and step outcomes, and uploads it even though the job failed.
The measured step requires all three preflight steps to have succeeded.

Before any measured child, the supervisor reruns #831’s pinned *harmless*
`harmless_cap_probe.py` on the actual target Actions image and archives its
stdout, stderr and JSON receipt. It requires hard=soft `RLIMIT_AS` at each
of the three inherited caps and an over-cap mapping rejection, with the
checkout head, Ubuntu image version and Python provenance bound to the
dispatch. This admits a later Ubuntu image only when that fresh control
passes; #831’s original preparation receipt remains independently pinned.
A probe failure consumes the one-shot label and is a pre-dispatch refusal.

The exact phase order is toy → n13 panel → n131, stopping at the first
failure. The inherited hard caps are unchanged:

| Phase | Outer wall | Hard address space | Other caps |
| --- | ---: | ---: | --- |
| toy | 60 s | 512 MiB | each CNF ≤1 MiB |
| n13 panel | 1,800 s | 2 GiB | each query CNF ≤8 MiB; solver 30 s/1 GiB; DRAT ≤128 MiB; checker compile/check 60 s each |
| n131 single edge | 90 s | 512 MiB | ≤500,000 DAG nodes; CNF ≤64 MiB |

The panel order is the nine frozen pairs, each SAT case immediately followed
by its wrong-output UNSAT case, then two invalid-point UNSAT cases. SAT
admission requires a *complete* clause-valid model lifted to the exact full
points and independent group sum. UNSAT admission requires a nonempty text
DRAT file, a bounded checker receipt with exit 0, and fresh checker replay on
the exact archived CNF. n131 replay rebuilds the frozen 131-bit relation,
checks all DAG metadata against the producer receipt, then streams the DIMACS and
compares every gate clause and the asserted output literal with the rebuilt
relation. Dimension-preserving gate and output mutations, plus the previously
accepted one-variable fake, are rejected by no-network controls. n131
admission remains only a single-edge export check under caps. It is not an
n131 PDP or attack claim.

The one cold run writes a new `evidence/run1/` archive after download. Its
`MANIFEST.json` hashes every raw file, including streams, per-case CNFs,
proofs, progress and phase results. `receipt.json` records exact release
provenance, commands, stop reasons, hard caps, wall, sampled RSS and phase
Python CPU; the latter excludes solver/checker subprocess CPU and must not be
presented as total process-tree CPU. Preserve both the #831 preparation/build
cost and v2 release/setup/phase costs separately. An Actions artifact alone
is not a durable record: commit the exact raw files and an outcome decision in
this or a linked PR, then require archive-only `ci_replay.py --evidence ...`
before accepting a pass. A failed or censored attempt is archived in full and
does not get rerun under the same release.

The release decision is therefore binary only for this *single-edge stage*:
`PASS_INDEPENDENT_REPLAY` after all three phases and fresh proof/model checks,
or `ARCHIVED_FAILURE_ONLY`/`ARCHIVED_PRE_DISPATCH_REFUSAL` with its exact
stop point. End-to-end ECC2K-130 factor-base, PDP and rho decisions stay with
the separately charged campaigns and canonical scoreboard.

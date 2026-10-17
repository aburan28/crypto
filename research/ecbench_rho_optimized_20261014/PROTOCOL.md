# Strong-rho implementation optimization protocol

Status: **preregistered; no measurements have run under this protocol**.

## Requested outcome and evidence contract

The requested outcome is to make the Pollard-rho arm used by ECBench's Koblitz
comparisons materially faster, without weakening the reference or changing the
problem it solves.  This round is an implementation-throughput optimization of
`rho.signed_frobenius_strong`, not a change to its generic-group boundary.

| Requirement | Acceptance evidence |
|---|---|
| Preserve the strong one-target reference | Same curve, public target, algorithm seed, 32 lanes, `dp_bits=4`, jump table, signed-Frobenius canonical representative, collision rule, and setup/online boundaries |
| Preserve deterministic work | For every matched run: identical recovered scalar, status, phases, charged GAE, operation counters, trail diagnostics, and automorphism count; all scalars verify as `[d]G=Q` |
| Establish a material native improvement | Candidate/baseline solve-wall ratio has an admitted paired 95% interval below 1, lies outside the optimized A/A interval, and the deterministic Callgrind solve instruction ratio is at most 0.85 on the frozen n=61 profile input |
| Preserve historical evidence | `ecbench verify --replay-all` is identical on the selected committed frontier sessions; any mismatch blocks the optimization |
| Keep the frontier honest | Recompute every affected same-target rho/IC comparison.  Faster rho can only reduce an IC-over-rho wall ratio; the leaderboard is regenerated when an admitted value changes and otherwise records why ranked rows stay fixed |

The 15% instruction target makes “materially faster” falsifiable.  A smaller
correct improvement may still be reported as engineering, but it does not meet
that target.  A count, answer, or replay change is inadmissible even when faster.

## Boundary and classification

For each curve, the derived generic floor is

```text
S_floor = sqrt(pi / (2 A)),  A = 2n,
S = cold group-addition equivalents / sqrt(r).
```

The measured reference is `rho.signed_frobenius_strong` at its immutable
registry defaults (`lanes=32`, `dp_bits=4`, `step_cap_factor=2000`).  The
unmodified implementation is retained under the new comparison-only method id
`rho.signed_frobenius_strong_legacy`; it must produce exactly the same counted
work.  `rho-strong-aa` repeats the optimized method and supplies the wall-time
noise interval.  Since the mathematical algorithm and `S` must not move, a
successful change is **engineering**, not an advance toward the generic floor.

## Frozen inputs

[`spec.json`](spec.json) contains six registered Koblitz curves already used by
the canonical medium-tier bound, eight independent public targets per curve,
one warm-up and five measured rounds.  Arms are alternated by ECBench in one
session.  No target, seed, curve, parameter, round count, timeout, accounting
field, or isolation threshold may be changed after a run starts.

The baseline source is commit `7badd44a3` on
`codex/ecbench-frontier-20261013`, including the Rust 1.98 lint follow-ups that
were not part of the already merged frontier PR.  Builds use Rust 1.98.1 and
the repository's default release profile.  Record the exact dependency
resolution, compiler, binary SHA-256, host capsule, CPU feature set, and dirty
state for both baseline and candidate.

The deterministic profile input is the first measured `rho-strong` execution
for the n=61 curve in the planned schedule.  Export it with `ecbench
profile-input`; retain the input JSON, its SHA-256, both binaries' child output,
Valgrind version, and every Callgrind part.  `callgrind.Ir` is a separate
implementation unit and never replaces `S`.

## Measurements

1. Build and archive the unmodified baseline binary before editing the hot
   path.  Record `ecbench host`, `rustc -Vv`, binary SHA-256, source commit and
   the plan for `spec.json`.
2. Profile the frozen n=61 child through the documented complete-solve
   Callgrind interval.  Attribute instructions by function before choosing an
   optimization.
3. Keep the literal pre-change canonicalization path for the `rho-legacy` arm;
   route the standing `rho-strong` id through the candidate path.  Both paths
   must return the same canonical point and coefficients for exhaustive edge
   cases and randomized curve states.
4. Run the session with the native harness at the highest level this managed
   host can earn, requesting L2.  Retain every warm-up, measured failure,
   timeout, isolation blocker and contended run.  Do not use `--allow-busy`.
5. Save `ecbench compare` results for legacy/optimized and optimized/A/A,
   `ecbench table`, and `ecbench resources`.  The primary runtime statistic is
   ECBench's paired wall ratio and bootstrap interval; the operation ratio must
   be exactly one.
6. Audit the new session with `verify --replay-all`.  Replay the committed
   strong-rho frontier sessions named in the result report with the candidate
   binary.  A different-host independent receipt remains required for any new
   `vs_rho` claim; this implementation benchmark alone is not such a claim.

No build, test, audit, or unrelated agent job may overlap a timed session.  The
host is a VM, so L3 is impossible; wall time is a result only if ECBench admits
L2 and the A/A gate passes.  Counts and Callgrind remain reportable if wall
admission fails.

## Success, stop, and rejection conditions

- **Correctness gate:** 100% of completed measured runs verify, all matched
  legacy/optimized count projections are identical, and all selected historical
  replays are identical.
- **Runtime gate:** at least five admitted pairs, candidate/baseline wall CI
  wholly below 1, outside the A/A interval, with no pooling across isolation
  levels.  Report per-curve rows; a pooled value is not a scaling claim.
- **Deterministic gate:** n=61 complete-solve Callgrind Ir ratio at most 0.85.
- **Stop/reject:** stop promotion on a wrong answer, changed trajectory/count,
  replay mismatch, missing setup or online phase, dirty measured binary,
  persistent timeout, or a candidate regression.  Preserve the attempt and its
  exact obstruction.

If the first candidate misses the instruction target, profile it and continue
only with changes supported by that profile.  Do not tune the frozen workloads,
seeds, or acceptance threshold to rescue the result.

## Required deliverables

- protocol, spec, plans, host/build manifests, frozen profile input and hashes;
- baseline/candidate source revisions and exact commands;
- sealed ECBench session, replay-all receipt, tables, comparisons, resources,
  unsuccessful attempts and historical replay receipts;
- a source-linked run report and result report, quantitative SVG graphs, PNG
  renders, inspected PDFs, and copies under `benchmarks/20261014-ecbench-rho-optimized/`;
- every affected frontier/scoreboard/leaderboard generated artifact in the same
  PR, or an evidence-backed statement that its canonical value did not change.

Public synthetic targets only.  This round makes no deployed-curve,
key-recovery, asymptotic, or ECC2K-130 claim.

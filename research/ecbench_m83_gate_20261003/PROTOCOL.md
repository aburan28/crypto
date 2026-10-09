# Protocol: the m = 83 gate curve in `ecbench` — step-rate diagnostic and the cost of the gate

Preregistered 2026-10-03, before any of the sessions below ran.

## Question

AGENTS.md §8a makes `E_0: y² + xy = x³ + 1` over `GF(2^83)` the confidence
gate every ECC2K-130-bound index-calculus improvement must pass, and §12
makes `ecbench` the harness such comparisons run in. Until this PR the
harness held `GF(2^m)` for `m ≤ 62` only. This protocol covers the first
thing a harness owes a curve it can now hold: **what one step of the
matched rho reference costs there, and therefore what a one-target solve
of the gate costs**, stated as an extrapolation from a measured step rate
and the walk's known floor, never as a solve.

It is a **stage diagnostic** (AGENTS.md §2, §6): no row below is a
verified solve, no `S` below is a method's `S`, and no speed claim
follows from it.

## Frozen inputs

| input | value |
|:--|:--|
| curve | `icv1-f2m83-tm6151469093347-debefd74` (`K_0 / GF(2^83)`), modulus `z^83 + z^45 + z² + z + 1` = `0x800000000200000000007` |
| subgroup | `r = 2417851639230796216685689` (prime), cofactor `4`, `#E = 9671406556923184866742756` |
| generator | `G = (0x477f77103dfad59850800, 0x2fa5e737d542c4e4fd5c3)`, the frozen generator of `research/ic_tool_program/conformance/v2/params/gate-m83-T001.json` |
| Frobenius eigenvalue | `λ = 254512724090651164922414` (`π(G) = [λ]G`, derived by `WideInstance::explicit` with the narrow curve's own routine and checked) |
| group | `koblitz_wide` (two-word field, affine arithmetic, the tuned walk ported statement for statement from `ic_boundary::rho_walk_with`) |
| comparison curve | `K_0 / GF(2^61)` built twice: narrow (`koblitz`) and wide (`koblitz_explicit` from the narrow plan's facts) |
| method | `rho.signed_frobenius_budget`, `steps = 20 000 000` walk operations per run (m = 83) and the same budget on the n = 61 pair |
| targets | 2 planted targets per curve (`target_seed` 83 and 61), 1 warm-up round, 3 measured rounds, measurement seed 2293761 |
| host | this Mac (Apple silicon, arm64, PMULL); `--cpus none`, so every run is L0: counts only, wall time descriptive |
| binary | the PR's `ecbench`, recorded in each session's `session.json` (`binary_sha256`, `git_commit`) |

Specs: [`specs/steprate-m83.json`](specs/steprate-m83.json),
[`specs/steprate-n61-narrow.json`](specs/steprate-n61-narrow.json),
[`specs/steprate-n61-wide.json`](specs/steprate-n61-wide.json).

## Boundary and reference

- The gate's floor for a one-target solve by a generic algorithm using
  the signed Frobenius (`A = 2n = 166`): `√(πr/2A) = √(πr/332) ≈
  1.513 × 10^11` walk steps, i.e. `S_floor = √(π/332) ≈ 0.0973`.
- The reference is the tuned walk itself; at every narrow size its mean
  `S` sits within a few percent of the floor
  (`research/ecbench_calibration_20261002`). No other method runs at
  m = 83 in this harness yet.

## Hypotheses

- **H1 (sameness, already a test):** on `K_0 / GF(2^61)` the wide
  construction records the same workload ids, the same counts and the
  same answers as the narrow one, record for record.  (`tests/ecbench.rs`
  checks this at n = 41; the n = 61 sessions here check it at the largest
  word-size degree, in a committed session.)
- **H2 (two-word overhead):** the wide walk's wall time per walk
  operation on n = 61 is between 1.3× and 4× the narrow walk's on this
  host.  Outside that band, something other than the extra word is
  being paid for and must be named before the number is used.
- **H3 (the gate's per-step cost):** the m = 83 wall time per walk
  operation is within 1.5× of the n = 61 wide figure (same word count,
  one more fold on the gate modulus).

## What is measured, and what is not

- Measured: walk operations, steps, distinguished points, canonicalisations
  and every tuned-walk counter per run; wall nanoseconds per run inside
  the measured process; host hardware counters where the host has them
  (none on this Mac).
- Derived: `ns / walk operation`; the ratio wide : narrow at n = 61.
- **Extrapolated, and marked so:** the gate's one-target cost,
  `E[steps] × ns/step` with `E[steps] ≈ 1.513 × 10^11 × (1 + δ)`, `δ`
  the tuned walk's measured overhead over the floor at narrow sizes
  (≈ 0.03–0.08 in the calibration sessions).  Reported in hours on this
  host, single-threaded, and labelled an extrapolation resting on H2–H3.

## Success and stop conditions

- The sessions complete; every run ends `exhausted` at or just past the
  budget (a solve inside 2 × 10^7 steps at m = 83 has probability
  ≈ 10^-8 per run and would be reported, not hidden); every record
  replays bit for bit (`ecbench verify --replay-all`).
- H1 holds exactly on the n = 61 pair (equal workload ids, counts,
  answers).
- H2 and H3 are read from the tables; a figure outside its band is
  reported as such with the suspected cause.
- Inadmissible: changing the budget, seeds or curve after seeing a
  number; quoting a wall figure from an L0 run as anything but
  descriptive; reporting the extrapolated gate cost as a measurement.

## The gate itself (pending)

A one-target solve of the gate is `≈ 1.5 × 10^11` walk operations.  At
the step cost this protocol measures it is hours per target on one core;
a panel of eight targets, the minimum the measurement skill asks of a
claim, is days.  That run is **pending**: it needs a dedicated host at
L2 or better for its wall figures (its counts are valid at L0), and the
strong single-target reference (`rho.signed_frobenius_strong`) is not yet
ported to two-word fields, so no IC-versus-rho claim at m = 83 can be
assembled by `ecbench claim` until it is.  Both are tracked in the PR's
closeout list.

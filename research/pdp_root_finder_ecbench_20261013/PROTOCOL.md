# Protocol: the root-finding decomposition oracles inside `ic.pipeline`, against the inherited-F4 oracle and the strong rho

Registered 2026-10-13, before any session of this directory ran.  Builds,
the unit tests of `ic_framework::rootfind` and `ecbench plan` on the specs
precede it; every measured run is produced after this file is committed to
the branch.  Companion: `research/pdp_root_finder_20261012/` (the oracles
measured alone, PR #1752) and `research/f6_ic_ecbench_ladder_20261005/`
(the inherited-F4 oracle on the same curves and bases, whose spec shape
this round copies).

## Question

Inside the repository's index-calculus pipeline, charged in ecbench's one
unit (`S` = group-addition equivalents / `√r`), what does an exact
root-finding decomposition oracle cost end to end — factor base, oracle
set-up, relations, linear algebra, verification — against (a) the
inherited-F4 Gröbner oracle on the same base and targets, (b) the harness's
own `mitm:m=3` table oracle, and (c) the strong single-target rho the IC
claim rules name as the reference; and does the comparison change with the
field degree?

## Boundary and reference

- Boundary: each curve's generic floor `√(π/4n)` in `S`, recorded per run.
- Reference arm: `rho.signed_frobenius_strong` (lanes 32, `dp_bits` 4,
  `step_cap_factor` 2000, as in the F6-IC ladder); `rho.signed_frobenius`
  as the operations-only baseline; `rho-strong-aa`, a second strong-rho arm,
  as the A/A wall-noise control.

## Arms (identical on every size unless stated)

| arm | role | method |
|:--|:--|:--|
| `rho-strong` | reference | `rho.signed_frobenius_strong` |
| `rho` | baseline | `rho.signed_frobenius` |
| `ic-f4` | baseline (`n ≤ 23`) | `ic.pipeline`, `koblitz-standard-subspace:dimension=d`, `pdp3-koblitz:m=3,engine=inherited-f4,degree=3,node_budget=8192` |
| `ic-mitm` | baseline | `ic.pipeline`, same base, `mitm:m=3,negation_folded=1` |
| `ic-rootfind` | candidate | `ic.pipeline`, same base, `rootfind-mitm:m=3` |
| `ic-rootfind-pairs` | candidate (`n ≤ 31`) | `ic.pipeline`, same base, `rootfind-pairs:m=3` |
| `rho-strong-aa` | control | `rho.signed_frobenius_strong` |

`max_trials = 200000`, `linalg = incremental-gauss` (the default), eight
public one-target workloads per curve, `warmup = 1`, `order = alternate`,
`isolation_required = L2`, 3,600 s timeout per run.

## Sizes

`E_0: y² + xy = x³ + 1` over `F_{2^n}`, the registered curves the F6-IC
ladder measured, with its base dimensions, plus `n = 37`:

| `n` | `d` | rounds | `ic-f4` | `ic-rootfind-pairs` |
|--:|--:|--:|:--|:--|
| 13 | 5 | 5 | yes | yes |
| 19 | 7 | 5 | yes | yes |
| 23 | 8 | 3 | yes | yes |
| 31 | 10 | 3 | no: inherited F4 hit the 3,600 s cap on its first target of this size in the F6-IC ladder | yes |
| 37 | 12 | 3 | no | no: about `A²` additions a trial at `A ≈ 2^11` abscissae, three minutes a run with nothing left to learn about a floor arm |

`K_1` (`a = 1`) was not used: the constructor has no usable prime-order
subgroup at `n = 13, 31, 37` and a 15-bit one at `n = 29`.

## Units and what is read

- `S` per arm and workload, the floor, and `S / S_rho` per (workload,
  round): the paired totals ratio `ecbench compare` reports, with its
  bootstrap interval.
- `ic-f4`'s `S` is a **lower bound** (the solver's word XORs are counted and
  not charged, as the F6-IC ladder read it); `ic-mitm`, `ic-rootfind` and
  `ic-rootfind-pairs` charge every group operation.  A ratio of a lower
  bound to a complete cost is read as a counted-work diagnostic, never as a
  speedup.
- Wall time only where the harness admits it (both runs at the spec's level,
  outside the A/A noise); otherwise descriptive.
- The replay audit (`verify --replay-all`) must be `OK` with every replay
  `identical`; a session that fails it reports no result.

## Predictions

P1. `ic-rootfind` and `ic-mitm` charge the same relation-phase additions on
    every workload (`2A` a trial, the same trial stream), so their `S`
    agree to the probe confirmations; their wall times differ.
P2. On every size with an `ic-f4` arm, `ic-rootfind`'s complete `S`
    exceeds `ic-f4`'s lower-bound `S` (the F4 arm's charged additions are
    the pipeline's, not the solver's), and `ic-rootfind`'s wall time per
    run is below `ic-f4`'s by at least 100× at `n = 23`.
P3. Every IC arm's `S` is above the strong rho's on every size (the F6-IC
    ladder: F4 1.1–2.1× rho at `n = 13…23` as a lower bound), and
    `S_IC / S_rho` does not fall with `n` for any oracle: no exponent signal.
P4. `ic-rootfind-pairs` charges about `A²` additions a trial and sits far
    above every other arm.

## Decision rule

- Audit `OK` on a size, else no result on that size.
- The per-size table reads the paired ratios; a pooled ratio across sizes
  is not a scaling claim.  A claim about growth needs the fit over the
  sizes that completed, with its interval, and at least four sizes.
- A `vs_rho` claim report is built for the first measured round of one
  workload per size (`ecbench claim build`); without an independent replay
  receipt from a host of another class it fails on exactly
  `independent_validation`, and that is reported as the claim's status.
- Class per AGENTS.md §3 from the ratio to the floor and to the reference:
  the expectation is **engineering** (a constant on the generic line, no
  exponent moved), with the honest end-to-end picture that index calculus
  with an exact oracle is above rho at every size measured here.

## Inadmissible

Changing the workloads, seeds, unit, base or accounting after seeing
results; pooling contended or below-level wall runs; quoting a lower-bound
ratio as a speedup; a claim built by any binary other than the one that
measured the session.

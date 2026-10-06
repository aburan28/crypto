# Bounds and frontiers

**What a method costs, as a record; which methods nobody has beaten, as a
page; how to beat one, as a procedure.**

`ecbench` measures (`run`), pairs (`compare`), tabulates (`table`) and checks
one IC-against-rho claim (`claim`). None of those records *the established
cost of a method* as an object another session can cite, re-derive or
challenge. The best known constant for a walk lived in prose and in the
registry's derived expectation; whether a change moved the exponent, the
constant, or only the accounting was left to the author. This directory adds
four sealed records and the commands that produce them:

| record | schema | id | what it says |
|:--|:--|:--|:--|
| **bound** | `ecbench.bound/v1` | `ECBND1h…` | the cost of one method on one domain: `S = gae / √r` at the declared exponent, `gae = C · r^α` fitted across sizes, both with 95 % intervals, per phase, with what it stored and what it left unpriced, read from named sessions and audit receipts |
| **frontier** | `ecbench.frontier/v1` | `ECFR1h…` | per domain, the bounds no other bound dominates; built from the records, never edited |
| **challenge** | `ecbench.challenge/v1` | `ECCH1h…` | what a candidate must do to replace an incumbent: domain, incumbent, curves, targets, rounds, acceptance rule, a nonce |
| **verdict** | `ecbench.verdict/v1` | `ECVD1h…` | one paired session judged against a challenge: `advances`, `trade`, `matches`, `regresses` or `inadmissible`, with the axes and the level it moved |

Everything below is in the harness's unit (`ecbench.gae`, group-addition
equivalents) on single-target workloads, under AGENTS.md's rules: the whole
method, cold, one target, operation counts as the metric, wall time as a
practicality note and never the basis of a bound.

## 1. What a bound is, and is not

A bound is a **measurement with its scope attached**: this method
configuration (`ECM1h…`), this problem, this curve family, these sizes, this
tier, this unit, these sessions, this many verified runs. It carries the
constant and the exponent *with intervals*, because rho's cost is a random
variable with a spread of the order of its mean, and a point without an
interval is a guess with a number on it.

A bound is **not** a lower bound on the problem, not a theorem, not a
statement about cryptographic sizes, and not a wall-clock figure. The generic
floor `√(π / 2A)` it is measured against is a model of collision search, not a
proof about implementations. **Toy stays toy**: a bound in the `toy` tier
(fields of at most 32 bits) says nothing about `medium` or `crypto`, and the
frontier never mixes tiers. This is `claim_tier` from
crypto-autoresearcher's `docs/claims-and-verification.md`, applied at the
record level rather than by convention.

## 2. Where a change can happen: four levels

Every claimed speed-up moves exactly one of these, and the protocol refuses
to let one masquerade as another.

| level | what moves | how it is read | unit |
|:--|:--|:--|:--|
| **exponent** `α` | the algorithmic class: `r^{1/2}` to `r^{1/3}` | the slope of `ln gae` against `ln r` over at least four sizes, with its interval disjoint from the incumbent's | any counted unit |
| **constant** `C` | walk design, automorphisms used, distinguished-point parameters, table shape | the paired ratio `Σ S_candidate / Σ S_incumbent` with its interval excluding 1, at the same `α` | `ecbench.gae` |
| **primitive weights** | formulas and coordinates: field multiplications, squarings and inversions per group operation | for the prime-field generic walks and tables (`rho.*`, `bsgs.*`, `kangaroo.vow` on `PrimeCurve`): the record's `field_ops` block, tallied inside the group law — every modular multiplication (`muls`, one by a small constant included), every squaring where the formula squares a value (`sqrs`), every modular inversion (`invs`) — read as the axes `field_muls`, `field_sqrs`, `field_invs` per `√r` (§5); **unknown**, absent and never zero, on binary and Koblitz curves and for the index-calculus pipeline | field operations |
| **machine** | cycles per primitive: micro-architecture, SIMD, memory layout | wall time, instructions retired, isolation level | seconds, Ir |

The example that motivates the table: *one fewer squaring in the point
addition formula*. It is a primitive-weight change. It cannot move `α`, and in
`ecbench.gae` — which charges a group addition as one unit whatever its field
cost — it cannot move `C` either. The `gae` figure of a bound is therefore
blind to it by design, and a report that turns it into "`r^0.9`" has confused
levels. It is real, and it is measured at its own level: on prime-field
curves the generic walks and tables count the field operations behind every
group operation they charge, the record carries them as `field_ops`, the
bound as the axes `field_muls`, `field_sqrs` and `field_invs` (§5), and a
verdict reads them paired. One fewer squaring in the addition formula is then
a move on `field_sqrs` at an unchanged `ops`, and the level it names is
`primitive`. Where nothing counts them — binary and Koblitz curves, and the
index-calculus pipeline, whose field arithmetic outside the group law (square
roots, Legendre symbols, oracle inversions, the elimination) is counted
natively and priced apart — the axes are unknown, not zero. The counts reach
the end-to-end figure only through composition (§7), as a *derived* number
until an end-to-end session in a field-operation unit confirms it.

The level a verdict names is a statement about operations: `exponent` when
both arms' fitted `α` intervals are disjoint over four or more sizes and the
candidate's is lower, `constant` otherwise when `ops` moved, and `primitive`
when `ops` did not move but a field-operation axis did. An advance on memory
alone names no level.

## 3. Records

Every record is pretty-printed JSON whose id is `prefix + 12 hex` of the
SHA-256 of **its own bytes with the id field empty**, the seal `Record::seal`
uses for runs. A changed byte is a different record; two sessions that fit the
same sessions write the same bytes; a tampered record fails `bound check`
before anything is read from it. Records are write-once: a correction is a new
record, and the old one stays. Timestamps are kept out of every record for the
same reason.

### 3.1 Domain

Bounds compare only inside one domain, and a domain's id (`ECDOM1h…`) hashes
every field of it:

```json
{"problem": "ecdlp.single_target", "family": "prime", "target_kind": "planted",
 "unit": "ecbench.gae", "tier": "toy",
 "envelope": {"targets": 1, "precomputation": "none", "threads": 1}}
```

`family` is `prime`, `koblitz` or `binary` (the automorphism group, and so the
floor, differs). `tier` is derived from the field size — `toy` ≤ 32 bits,
`medium` ≤ 96, `crypto` above — and a fit whose sizes span tiers is refused
(`--tier` keeps one). `envelope` is the resource shape every ecbench workload
has today; a multi-target or precomputation method is a different domain, not
a better entry in this one.

### 3.2 Bound

| field | meaning |
|:--|:--|
| `method` | the configuration: registry id, resolved parameters, `method_id` (`ECM1h…`). The code's hash is in the provenance, not the identity: a registered method never changes what it counts (CI replays every committed session), so the configuration names the counts. |
| `level` | `exponent` when four or more sizes verified (a scaling claim), `constant` otherwise |
| `sizes[]` | one row per curve, sorted by `r`: slug, `log2_r`, field bits, `A`, the floor, workloads, runs, verified runs, mean `gae`, mean `S` with its interval, ratio to the floor, the registry's derived `S` and its ratio to the floor, memory and unpriced work per `√r`, isolation levels earned |
| `fit` | `gae = C · r^α`: `alpha` and `alpha_ci95`, `log2_c`, `r_squared`, the declared `α` (`0.5` for every method the registry derives an expectation for), whether the interval contains it, `scaling_claim` |
| `constant` | mean `S` and mean ratio to the floor with intervals; the derived ratio to the floor when it is one number across sizes (`1.0` for `rho.negation`; `1.5 / √(π/4) = 1.69` for `bsgs.textbook`) |
| `dimensions` | the axes (§5): `ops`, `memory`, `uncharged`, each with `known`, `value`, `ci95`, its statistic and source; and `field_muls`, `field_sqrs`, `field_invs` (mean per `√r`, the same two-stage interval) **only when every verified run of the arm carries `field_ops`** — otherwise the keys are absent, not unknown, so a bound fitted from records written before the block existed is byte for byte the bound it was |
| `stages[]` | per phase name: mean `gae`, share of the total, sizes with work, its own fitted `α` with interval when it has work on four or more sizes |
| `provenance` | sessions (directory, session id, spec id, status, arm, role, binary hash, commit, environment class, SHA-256 of `records.jsonl`), audit receipts (path, hash, `ok`, replays, replays reproduced), record counts by status |
| `admissibility` | `admissible` or `inadmissible` (a measured run did not verify, or a receipt reports problems), `bounded` (some work was counted and not priced: every total is a floor), the unpriced counters, determinism, reasons |
| `fit_options` | the filters and bootstrap settings, so `bound check` fits it again |
| `improves_on`, `verdict_id` | set by the verdict that produced the bound; empty on a bound fitted from a session directly |

Every statistic carries a two-stage bootstrap interval — sizes resampled,
then runs within each — the interval `compare` and `table` already use, with
2 000 resamples from a fixed seed so a bound re-derives bit for bit.

The run records a bound reads (`ecbench.record/v1`) carry one optional block
of their own for the primitive level: `field_ops: {muls, sqrs, invs}`, the
modular multiplications, squarings and inversions behind every addition and
doubling of the solve, summed over its phases. It is written only when the
method's group counted it (prime-field curves through the generic walks and
tables), so an older record is unchanged and its seal still checks; it lives
outside `counters` and `phases`, which the replay of a committed record
compares exactly; and `verify` compares it on replay only when both the
record and the replay carry it. What is counted: every modular multiplication
the group law performs is a `mul` (the doubling's multiplications by the
constants 3, 2 and 2 included, since the code performs them as modular
multiplications), every squaring where the formula squares a value is a
`sqr` (so `sqrs` is a genuine subset, not a reading of equal operands), every
modular inversion is an `inv`; additions, subtractions, negations and
comparisons are not counted. In the code as written an affine addition is
`2M + 1S + 1I` and a doubling `5M + 2S + 1I` (`PrimeCurve::FIELD_OPS_PER_ADD`,
`FIELD_OPS_PER_DOUBLE`); the special cases `∞ + P` and `P + (−P)` do no field
arithmetic. Binary and Koblitz curves and the index-calculus pipeline write no
block: unknown, never zero.

### 3.3 Frontier

Built by `frontier build` from a directory of bounds. Per domain: every
admissible entry with its axes, `is_frontier`, `dominated_by`, the reasons,
the ops and memory leaders. Inadmissible bounds are listed and never on a
frontier. The page `FRONTIER.md` and `frontier.json` are generated and
committed; CI rebuilds both and fails when they are stale.

### 3.4 Challenge

| field | meaning |
|:--|:--|
| `domain`, `domain_id` | the domain the challenge is in |
| `incumbent` | the method to beat, as a method spec, and its committed `bound_id` when it holds one |
| `workloads` | the curves (at least `acceptance.min_sizes`, all in the domain's family and tier — checked by building them), targets per curve, target kind, and a `nonce` |
| `measurement` | rounds, warm-up, required isolation level, timeout |
| `acceptance` | `axes` dominance reads (default `ops`, `memory`; `uncharged`, `field_muls`, `field_sqrs`, `field_invs` opt-in); `min_sizes` (4); `min_runs_per_size` (8); `require_audit`; `require_replay_all`; `uncharged_tolerance` (0.05) |

The spec a candidate runs is a function of the challenge, the candidate's
method spec and an **epoch**: `challenge spec --epoch N` derives the target
seed and the algorithm seed from the nonce and `N`, and writes three arms —
`incumbent` (reference), `candidate`, `incumbent-aa` (an A/A control of the
incumbent) — interleaved on the same workloads. Anyone can rebuild the spec and
check its `ECS1h…` id; nobody can tune to targets that do not exist until an
epoch is named.

### 3.5 Verdict

| field | meaning |
|:--|:--|
| `session` | directory, session id, its spec id, the spec id the challenge yields for this epoch and candidate, whether they match, status, binary hash, commit, levels earned |
| `audit` | the audit run by the verdict itself: `ok`, problems, records, verified, replays and how many reproduced, whether every run was replayed, and the SHA-256 of every session file |
| `incumbent`, `candidate`, `control` | arm names and method ids; the control's A/A ratio and interval (`1.000 [1.000, 1.000]` for a deterministic method under shared seeds) |
| `axes[]` | per axis: both means, `Σ candidate / Σ incumbent` over matched `(workload, round)` pairs, its interval, pairs, `better` / `worse` / `indistinguishable` / `unknown`, and whether the axis decides; `ops`, `memory`, `uncharged` and the three field axes are always reported, the field axes `unknown` unless every measured run of both arms carries `field_ops` |
| `per_curve[]` | the ratio per size, the rows a scaling claim reads |
| `fits` | both arms' `α` with intervals, and `exponent_moved` |
| `stages[]` | per phase: each arm's share and the paired ratio — which sub-algorithm moved (§7) |
| `accounting` | unpriced counters on each side, `bounded`, the unpriced-work ratio and `uncharged_shift` |
| `incumbent_drift` | the incumbent's recorded ops figure against what it measured in this session, when `--bounds` is given: a disagreement is a reason for inadmissibility, because the recorded bound and the session cannot both be right |
| `outcome`, `advances_on`, `regresses_on`, `level_moved`, `reasons`, `statement` | the decision and why; `level_moved` is `exponent` or `constant` when the advance includes `ops`, `primitive` when it does not and includes a field axis, `null` otherwise |

## 4. The fit

Points are every verified measured run of the arm, `(ln r, ln gae)`, one
stratum per size. Least squares gives `α` and `ln C`; `R²` says how much of the
spread the law explains. The interval is the two-stage bootstrap of the slope.
`declared_alpha` is `1/2` wherever the registry derives an expected `S`
(`methods::expected_s`), and `alpha_agrees_with_declared` says whether the
interval contains it.

Three things to know when reading `α`:

- **Fewer than four sizes is not a scaling claim** (AGENTS.md §5). The bound
  then has `level: constant`, and `α` is printed as description.
- **Fixed costs pull `α` below the law at small sizes.** A method whose
  set-up (jump table, lanes, factor base) is a large share of the total at
  `2^16` fits a flat line there. The stage rows carry the set-up's share and
  the search phase's own `α`; read those before reading the headline.
- **Rho's variance is the law's.** With 24 runs per size on four sizes,
  `rho.negation`'s `α` interval is `[0.27, 0.61]`. That is what the data
  supports, and the record says so instead of rounding to `0.5`.

The constant at the declared exponent is `S = gae / √r` averaged over verified
runs, and its ratio to the floor `√(π / 2A)` of each curve. On prime curves
`A = 2` and the two are proportional; on Koblitz curves `A = 2n` varies with
the size, so the ratio to the floor is the figure that compares across sizes
and across families. It is the `ops` axis.

## 5. Axes and dominance

| axis | statistic | source | default |
|:--|:--|:--|:--|
| `ops` | mean `S / √(π / 2A)` | `cost.total_gae`, `r`, `A` | decides |
| `memory` | mean table entries per `√r` | the first of `inserts_uncharged`, `table_inserts_uncharged`, `table_entries`, `distinguished_points` the record carries | decides |
| `uncharged` | mean `Σ *_uncharged` counters per `√r` | every counter the unit counts and does not price | reported; decides when a challenge names it |
| `field_muls` | mean modular multiplications per `√r` | `field_ops.muls` | on a bound only when every verified run carries the block; on every verdict; decides when a challenge names it |
| `field_sqrs` | mean modular squarings per `√r` | `field_ops.sqrs` | as `field_muls` |
| `field_invs` | mean modular inversions per `√r` | `field_ops.invs` | as `field_muls` |

**Unknown is not zero.** A method that reports no table counter has `memory:
known = false`; the axis is left out of every comparison involving it and the
entry says so. The index-calculus pipeline is in that position today.

The three field axes are the primitive level (§2). A bound carries them only
when every verified run of its arm carries `field_ops` — today the generic
walks and tables on prime-field curves — and a bound without them has no such
key at all, so the committed records fitted before the axes existed are
unchanged. On a frontier an axis one entry lacks is left out of comparisons
with it, as `memory` is for index calculus, and the page shows the columns in
a domain only when some entry there has them. In a verdict they are reported
always and `unknown` unless every measured run of both arms carries the block.
They are not a unit: `ops` still decides in `ecbench.gae`, and a candidate
that is clearly better on `field_sqrs` and indistinguishable on `ops` has
moved the primitive level, not the constant.

On a frontier, bounds were measured apart and cannot be paired, so dominance
is conservative. `A` dominates `B` when, on every axis both know, `A`'s point
estimate is no larger and `A`'s interval does not lie wholly above `B`'s, and
on at least one axis `A`'s interval lies wholly below `B`'s. Ties stand: a
domain can have several frontier entries, and the one-sided answer is a
challenge. The point condition keeps a wide interval from pushing a precise,
better entry off the frontier through a secondary axis.

In a verdict the arms *are* paired, so each axis reads the interval of
`Σ candidate / Σ incumbent` against 1: `better` below, `worse` above,
`indistinguishable` across. The outcome is Pareto over the deciding axes:

| outcome | meaning |
|:--|:--|
| `advances` | better on at least one deciding axis, worse on none |
| `trade` | better on some, worse on others (BSGS against rho: fewer operations, a `√r` table) |
| `matches` | better on none, worse on none |
| `regresses` | worse on some, better on none |
| `inadmissible` | the session is not the challenge's spec for the epoch; the audit found a problem or did not replay every run when required; too few sizes or runs; a measured run did not verify; the incumbent's recorded bound disagrees with its measurement here. The paired ratio is still reported, for information. |

## 6. Challenge to verdict

```sh
# 1. The challenge, from a draft (domain, incumbent, curves, measurement,
#    acceptance); `seal` builds the curves, checks family and tier, hashes.
ecbench challenge seal --draft draft.json --out docs/bounds/challenges/prime-toy-ops.json

# 2. The spec for your candidate in a fresh epoch (any integer nobody has used
#    for this challenge; the verdict records it).
ecbench challenge spec --challenge docs/bounds/challenges/prime-toy-ops.json \
    --candidate '{"id":"kangaroo.vow2","params":{"cap_multiple":"64"}}' --epoch 3 --out spec.json

# 3. Run it as any session is run (docs/ecbench/README.md; L2 for wall time,
#    irrelevant to the counts).
ecbench run --wait --spec spec.json --out research/<topic>/sessions/epoch3

# 4. The verdict: audits with every run replayed, pairs, fits, decides; writes
#    the candidate's bound on an advance or trade.
ecbench challenge verdict --challenge docs/bounds/challenges/prime-toy-ops.json \
    --dir research/<topic>/sessions/epoch3 --epoch 3 --replay-all \
    --bounds docs/bounds/records --out research/<topic>/verdict.json \
    --bound-out docs/bounds/records/<written by you>.json --audit-out research/<topic>/audit.json --exit-code

# 5. Rebuild the frontier; commit the session, the verdict, the receipt, the
#    bound and the page in one PR.
ecbench frontier build --bounds docs/bounds/records --out docs/bounds/frontier.json --markdown docs/bounds/FRONTIER.md
```

A candidate is a registered method (`ecbench-extend`): a changed algorithm is
a new registry id, the old id keeps its evidence, and the challenge's spec
names both. A candidate binary against a baseline binary is what two sessions
of one spec and `compare --b-dir` are for; a verdict reads one session, so a
code change goes through the registry.

The verdict is **re-derivable**: run it again on the same session and you get
the same `ECVD1h…`. It embeds the SHA-256 of every session file, not the
receipt's own hash (which carries a timestamp).

## 7. Sub-algorithms: stages and composition

A method is phases — `setup`, `search`, `internal_verification` for the tuned
walks; `factor_base`, `oracle_setup`, `relations`, `linear_algebra`, `verify`
for index calculus — and every record already charges each phase separately.
A bound carries each phase's share of the total and, with work on four or more
sizes, the phase's own `α`; a verdict carries each phase's paired ratio. So a
reader sees *where* a change landed: a new jump-table construction shows in
`setup`; a better walk in `search`; a candidate whose `relations` share fell
while `linear_algebra` rose traded one sub-algorithm's cost for another's.

Composition across levels is a sum over phases of counts times weights:

```text
cost_u(method, r) = Σ_phase Σ_op  count(op, phase, r) · w_u(op)
```

`w_gae` is 1 for an addition or a doubling and 0 for everything else (the
unit's charging rule, `docs/ecbench/README.md` §7). A field-operation unit
would weight an addition by its multiplications, squarings and inversions on
the curve family at hand, and a time unit by pinned nanoseconds per native
operation (`docs/ic/calibration.json`, the `ICBCAL1h…` discipline of
`aburan28/cryptanalysis`). On prime-field curves the weights are measured
rather than assumed: every phase's `adds` and `doubles` and the solve's
`field_ops` (§3.2) are the raw material of `w_field`, since dividing the
solve's multiplications, squarings and inversions by its additions and
doublings gives the weight per group operation the code actually paid —
`2M + 1S + 1I` per affine addition and `5M + 2S + 1I` per doubling as written
today — and a candidate formula shows up as a different quotient on the same
group-operation counts. **A composed figure is derived, never measured, and
never enters a frontier**: it predicts what an end-to-end session in that
unit should find, and the session is what moves the frontier in that unit.
This is where "one fewer squaring" lives — as a weight change, now visible on
the `field_sqrs` axis of the records that count it, whose end-to-end
consequence in any other unit is a prediction until measured — and it is why
the four levels are kept apart.

## 8. Rules

1. **Records are write-once and content-addressed.** A changed bound is a new
   bound. Nothing overwrites; `frontier.json` and `FRONTIER.md` are generated
   from the records and checked in CI.
2. **A bound is re-derivable.** `bound check` fits it again from the sessions
   and receipts it names and requires the same id; CI does this for every
   committed record. A bound whose sessions are missing, edited, or re-count
   fails.
3. **One domain per bound, one tier per domain.** Mixed families or tiers are
   refused at fit time. A `toy` frontier is a statement about toy sizes.
4. **Wall time never enters.** No axis, no fit, no acceptance rule reads a
   clock. The isolation level is recorded for the reader of wall figures
   elsewhere.
5. **Unknown is not zero, and unpriced is not free.** A missing counter leaves
   an axis unknown; counted-but-unpriced work is reported on every record and
   verdict (`bounded`, `uncharged`), and a challenge may make it decide.
6. **A frontier is conservative; a verdict is paired.** Ties stand on the
   frontier. Only a paired session on the challenge's frozen workloads, audited
   with every run replayed, moves it.
7. **Levels are named, never inferred from the headline.** `exponent` needs
   disjoint `α` intervals over four or more sizes. Everything else that moves
   `ops` at fixed `α` is `constant`. `primitive` needs a field-operation axis
   to move while `ops` does not, and exists only where the axes are known
   (prime-field generic walks and tables). Machine changes are not measured
   here and are not claimed.
8. **Everything cites.** A bound names its sessions by directory, session id
   and `records.jsonl` hash, and its receipts by hash; a verdict names the
   session's files by hash and the challenge by id; a new bound names the
   verdict and the bound it improves on.

## 9. Reading the seed frontier

`records/` holds the bounds fitted from two committed sessions
(`research/ecbench_all_candidates_20261003/sessions/{prime,koblitz}`: four
sizes each, eight targets, three rounds, audited with 12 replays). What they
say, read honestly:

- **Prime, toy.** `bsgs.negation` leads on operations (`1.154 ×` the floor,
  `[1.048, 1.254]`) with a `0.5 √r` table; `rho.negation` is the memory-light
  walk at `1.466 ×` the floor, `[1.278, 1.653]`, storing `0.034 √r` points.
  Both are on the frontier: the time–memory trade, as `docs/ecbench/README.md`
  §7 says, and not a finding. So is the frozen reference walk
  (`rho.frozen_reference`, `4.066 ×`), through memory alone: it stores the
  fewest distinguished points of any entry, and the rule does not let a worse
  operations figure knock a better memory figure off. `rho.plain` (`1.990 ×`)
  is dominated by `rho.negation`: the negation map's `√2`, measured, with the
  canonicalisation it costs counted under `uncharged` and shown beside it.
- **Koblitz, toy.** Three entries stand: `bsgs.negation` on operations
  (`5.877 ×` the floor `√(π / 4n)`), `rho.signed_frobenius` on memory
  (`0.016 √r`), and `rho.negation` between them. The strong lockstep rho is
  dominated at these sizes because its set-up is 72 % of its total and its `α`
  fits at `0.13` — a fixed cost at `2^16` to `2^21`, not a law. Its `search`
  phase alone is the number to read, and the record carries it. The signed
  Frobenius walk's own operations interval, `[2.9, 19.9]`, is as wide as the
  four curves' cofactors are different; it holds its place on memory, not on
  a precise operations figure.
- **Index calculus** is inadmissible on the Koblitz session (half its runs did
  not verify at the two largest sizes) and far from the floor on the prime
  session (`22.7 ×` and `1258 ×`). Both are listed; neither is on a frontier.
  Its memory is unknown (no table counter), which is a gap in the pipeline's
  reporting, not a strength.
- **Every `α` interval contains 0.5** except where set-up dominates. Nothing
  in these records is a scaling result; they are constants with their scope.

## 10. Commands

```text
ecbench bound fit --dir SESSION... --arm ARM [--tier T] [--curve SLUG...] [--audit RECEIPT...]
                  [--label L] [--notes N] [--root REPO] --out BOUND.json
ecbench bound check --record BOUND.json... [--root REPO]
ecbench frontier build --bounds DIR|FILE... [--axes ops,memory] --out F.json --markdown F.md [--check]
ecbench challenge seal --draft D.json --out C.json
ecbench challenge check --file C.json...
ecbench challenge spec --challenge C.json --candidate '{"id":...}' --epoch N --out spec.json
ecbench challenge verdict --challenge C.json --dir SESSION --epoch N [--replay-all] [--bounds DIR]
                          [--root REPO] --out V.json [--bound-out B.json] [--audit-out A.json] [--exit-code]
```

`--axes` and a challenge's `acceptance.axes` may name `ops`, `memory`,
`uncharged`, `field_muls`, `field_sqrs` and `field_invs`.

`refit.sh` regenerates every seed record and the frontier from the committed
sessions; CI runs `bound check` on every record, `challenge check` on every
challenge and `frontier build --check` on the page.

## 11. Layout

```text
docs/bounds/
  README.md          this protocol
  FRONTIER.md        generated page, checked in CI
  frontier.json      its machine-readable twin
  records/*.json     bound records, write-once, named by id
  challenges/*.json  standing challenges, one per domain and axis leader
  refit.sh           regenerate records and frontier from committed sessions
```

## 12. Beyond this repository

- **crypto-autoresearcher** reads bounds into its ledger through an
  evidence record's `measured_bound` block, reviewed by the Coordinator like
  any evidence; a frontier move there is a decision, not a push
  (`docs/bounds-and-frontiers.md` in that repository). The tier rules are the
  same rules.
- **cairn** holds the same frontier as a ratchet objective per domain, with a
  pinned evaluator that rescores a bound record from its own per-size counts
  and refuses a wall-clock basis; the network's `frontier_status` and this
  page then agree or one of them is wrong.
- **aburan28/cryptanalysis** fits the same law (`ca_bench complexity`) in its
  own units; a bound record in this schema from its receipts is the way to
  put the two libraries on one page without pretending the units agree.

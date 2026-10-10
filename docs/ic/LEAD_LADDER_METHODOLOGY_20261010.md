# Lead ladder for ECC2K-130 index calculus: baseline and measured gains

**Date:** 2026-10-10
**Conductor task:** T-248
**Status:** methodology (binding for the lead rows under
`research/ic_leads_20261010/`), lead registry, and the first baseline run.
Nothing here is a claim about the challenge curve. Public synthetic
known-answer fixtures only.
**Companions:** `RESEARCH_ECC2K130_IC_FEASIBILITY.md` (where the ladder stands
at n = 131), `PLAN_IC_ACCOUNTING_FIXES_20261007.md` (what may and may not be
called a speedup), `measurement/README.md` (ICMS, the comparability standard),
`RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md` and
`RESEARCH_ISOGENY_CONDUCTOR_GAP_20261007.md` (the leads this ladder measures).

## 0. The problem this solves

The repository has many index-calculus (IC) producers, several units called
`S`, and leads recorded as hypotheses in three research notes. A lead chased
in one producer cannot be read against a baseline taken in another, and a
wall-clock figure from this host (load average 90 to 150 during this session)
cannot be read against anything. This document fixes **one baseline
configuration per rung**, **one metric**, and **one protocol** under which
every lead is a declared single-field change against that baseline, so that
"lead X bought Y" has a meaning that survives the host and the author.

It does not replace the one-target `vs_rho` ledger
(`BOUNDARY_TARGETS.md`), which answers a different question (does IC beat rho
online on one frozen point). The ladder here answers: **by how much does each
structural idea move each phase of the method, in operations, at fixed
instance**.

## 1. Metric

The primary metric is the **operation count in `crypto.S.gae_pinned`**
(group-addition equivalents with the ratios pinned in `calibration.json`), as
`ic bench` reports it, per phase and in total, normalised by `√r` to give `S`
(`FRAMEWORK.md` §2). Reasons:

- Its group-operation part (additions, doublings, scalar multiplications per
  phase) is a deterministic count for every component used here (table
  oracles and Gauss over `Z/rZ`); no solver is priced from wall time in the
  baseline. **Caveat found in the first run:** `ic bench` reports
  `calibration_pins.pinned = []` on the `K_0` instances, so the native
  counters that are not group operations (Artin–Schreier solves in the base
  build, `row_ops` in the linear algebra, lookups) are converted to GAE at
  ratios *measured on this host* at the time of the run. Two identical runs
  gave LA GAE of 0.22 and 8.7 for the same 225 row operations. The lead
  tables therefore carry the raw counters (`fb_as_solves`, `setup_adds`,
  `decomp_adds`, `la_row_ops`, …) beside the GAE, and a ratio is read from
  the counters; the GAE column is comparable only once the instance's ratios
  are pinned in `calibration.json`.
- It survives the host. This session's host cannot earn ICMS isolation level
  L1, so **no wall-clock figure from this session is a result**; wall time
  is carried as a practicality note only.
- Every phase is inside it: factor base, table setup, failed decompositions,
  linear algebra, verification.

Beside the total, every row carries the **structural metrics** the leads claim
to move, because a lead that moves the total without moving the metric it
claims has been misattributed:

| metric | source field | which leads claim it |
|---|---|---|
| `columns`, `points_per_column` | factor base | folding (L1), invariant base (L4) |
| `rows` (relations needed), `rank` | linear algebra | folding, invariant base |
| `targets_tried / relations_found` (trials per relation) | decomposition | τ-slots (L5, L6), residual matching (L3) |
| `canonicalisations`, `frobenius_maps`, `lookups` (native, uncharged) | oracle counters | canonicalisation (L2) |
| `pair_table_entries`, `pair_table_representatives` | `mitm-frobenius-counted` setup | folding, canonicalisation |
| `n_vars`, `solving_degree`, Macaulay size, conflicts | algebraic oracles | Galois-orbit / equivariant solving (L7) |

Two derived quantities are reported on every row and are the only ones that
may be read across rungs:

- **`S` and `× reference`**: total operations over `√r`, and the ratio to the
  counted matched rho on the same instance and planted targets
  (`rho.measured_matched`, 16 runs). On a Koblitz curve the reference walk
  folds by negation and Frobenius (`A = 2n`); the floor beside it is
  `√(π r / (4n))` steps.
- **`B·D²/N`** (catalogue A1-4, the only certificate of non-genericity): `B`
  the signed point count of the base, `D` the online cost per target in the
  same unit, `N = r`. Every row in the current ledger has `B·D² ≫ N`; a row
  below 1 would be the first non-generic signal and must be replicated before
  it is reported.

## 2. Baseline configuration

One frozen configuration per rung, on the ECC2K-130 family member `K_0`
(`a = 0`, `y² + xy = x³ + 1`), the same curve shape as the challenge:

| field | value | why |
|---|---|---|
| curve | `K_0` over `F_{2^n}` | the challenge curve's family |
| factor base | `koblitz-orbit:divisor=1` (the smallest non-trivial Frobenius-invariant subspace, trace-zero), signed Frobenius orbits folded | the structure the leads are about |
| oracle | `mitm-frobenius-counted`, `m = 3` | the production shape of the landed rungs (pair table + one addition), with the native setup counters on |
| targets | `walk` (r-adding walk, one addition each) | the engineering ledger's `_walk` row |
| linear algebra | `incremental-gauss` over `Z/rZ`, stop at first log | the framework default |
| reference | counted matched rho, 16 runs, `A = 2n` | `rho.measured_matched` |
| repeats | 8 fresh planted targets per configuration, seed `20261010` | the `ic bench` ladder's noise floor; yield claims need ≥ 30 |
| unit | `crypto.S.gae_pinned` | §1 |

**Rungs.** Invariant subspaces exist only when `ord_n(2) < n − 1`; the table
below is the whole menu for `n ≤ 131`, and it is why the orbit-union base
(`compact-orbit-scan`) rather than a subspace is the only option at the
landed rungs `n = 53, 61, 83` and at the challenge `n = 131`:

| n | `ord_n(2)` | invariant subspace dimensions | prime `r` on | ladder role |
|---:|---:|---|---|---|
| 13 | 12 | 12 | `K_0` | baseline rung B0 (smoke; `2^12` abscissae) |
| 17 | 8 | 8, 16 | `K_1` only | baseline rung B1 on `K_1` (`K_0` has `r = 137·239`) |
| 19 | 18 | 18 | both | pair table too large (`2^18` abscissae) |
| 23 | 11 | 11, 22 | both | baseline rung B2 (`2^11` abscissae) |
| 29, 31, 37, 43, 47 | 28, 5, 36, 14, 23 | various | neither | no prime-order Koblitz instance at either `a` |
| 41 | 20 | 20, 40 | `K_0` | subspace has `2^20` abscissae: no pair table fits, so this rung belongs to the orbit-union arm |
| 53, 59, 61, 67, 83, 101, 107 | n − 1 | none | — | **orbit-union only** (compact-orbit pipeline) |
| 73 | 9 | 9, 18, 27, … | — | the one wide rung with small subspaces |
| 89 | 11 | 11, 22, 33, … | — | same |
| 127 | 7 | 7, 14, 21, … | — | same |
| **131** | **130** | **none** | `K_0` | **orbit-union only** |

The subspace arm therefore has three rungs, `n = 13, 17, 23`, and the
honest reading of the table is that the Frobenius-stable subspace is a toy
object: it coexists with a prime-order Koblitz subgroup and a buildable pair
table only below 24 bits. Everything the leads say about `n = 131` has to be
carried by the orbit-union arm.

The `ic bench` ladder (B1–B4) measures the leads where the production
pipeline's components are all counted. The compact-orbit pipeline
(`koblitz_rank_fixture`, `n = 41, 53, 61, 71, 73, 83`) is the second arm; its
baseline is the landed `K = 600` configuration of each rung
(`BOUNDARY_TARGETS.md`), and leads that only make sense there (A1 multiplier
closure, implicit S3 index, 5-sums) are measured as probes per relation and
states per column against that configuration, with the same declared-change
rule.

## 3. Protocol for measuring a lead

1. **Declare the change.** A lead run differs from the baseline in exactly
   one spec field (a plug-in name, one parameter, or one solver option). Two
   changes are two rows. Instance, seed, target law, unit, reference, stop
   rule and thread count may not differ (ICMS "forbidden" fields).
2. **Pre-register the prediction** in the lead's row of §4 before the run:
   which structural metric moves, by what factor, and which must not move.
3. **Run interleaved**: baseline, lead, baseline (A/B/A) in one `ic bench`
   sweep so both arms see the same planted targets (same seed, same walk).
   The second baseline arm is the noise floor for the row.
4. **Report** the ratio lead/baseline per phase and in total, the structural
   metrics of §1, and the verdict mix (found / refuted / gave up / wrong).
   A row whose answer did not verify is a failure row, never dropped.
5. **Classify** the change under `AGENTS.md` §3: *advance* (the exponent or
   the relation shape changed), *engineering* (constant factor on a phase),
   *relabelling* (the same work under a different column map), or
   *accounting* (a charged phase moved). The folding rows are the standing
   example: the `2n` is credited to a generic algorithm by the floor, so it
   is engineering, not an advance.
6. **Yield claims** (trials per relation, hit rate) need ≥ 30 targets and a
   binomial interval; **solver claims** need the per-call verdict and cost
   distribution, not the total; **wall-clock claims** need an ICMS session at
   the spec's isolation level, which this host cannot provide.

A lead "lands" on this ladder when its row has a measured ratio at two rungs
with the predicted structural metric moving and the others not, and the
classification is written down. Landing here does not promote any ledger
row; it is the input to one.

## 4. Lead registry

Each lead as the user named it, what it means in this codebase, where it
stands, the toggle that measures it, the pre-registered prediction, and the
falsifier.

### L1. Signed Frobenius folding

**Meaning.** Each signed Frobenius orbit `{±π^k P}` (`2n` points) is one
unknown; a relation among points becomes a row over orbit columns with
coefficients `±λ^k`. **Status:** implemented (`koblitz-orbit`,
`mitm-frobenius`, `_frobfold`); the `2n` is already in every landed rung.
**Toggle:** `koblitz-orbit:divisor=1` vs `koblitz-orbit:divisor=1,no_fold=1`
(base), `mitm-frobenius` vs `mitm` (table). **Prediction:** `columns ÷ 2n`,
`rows ÷ 2n`, LA operations `÷ (2n)²`; decomposition operations per target
unchanged within noise; table entries `÷ 2n`. **Falsifier:** a realised fold
below `n` (cofactor-projection merges can take it below `2n` but not below
`n`), or a decomposition-phase change above 20 %. **Class:** engineering (the
floor already credits `A = 2n`).

### L2. Canonicalisation

**Meaning.** The orbit representative of an abscissa under `x ↦ x²` and
negation is found by a normal-basis canonical form (minimum cyclic rotation of
the ONB coordinates); every probe pays one, and the table build pays one per
pair. **Status:** implemented and counted (`canonicalisations`, uncharged by
the unit). **Toggle:** none exists yet; the lead is a cheaper canonical key.
Two candidates: (a) an **orbit-invariant hash key** (a Frobenius-invariant
function of `x`, e.g. the coefficients of the minimal polynomial of `x` over
`F_2`, or a trace tower) used as the table key, with the rotation recovered
only on a hit; (b) a **partial canonical form** (compare on the top word of
the rotation only, full rotation on tie). **Prediction:** the count of
canonicalisations is unchanged (one per probe), the **cost per
canonicalisation** drops (ns/canonicalisation, reported by the counted
oracle's native timer), and nothing else moves. **Falsifier:** any change in
hit rate or relation count (the key is then not a complete invariant), or a
cost drop below 1.5× (then the rotation was not the bottleneck). **Class:**
engineering. **Honest note:** `KIC_TARGET_THREADS` showed the online stage is
probe-rate bound (≈ 1.1 M probes/s single core, dominated by the S3 solve,
not the canonicalisation), so this lead's ceiling is the canonicalisation's
share of a probe, which the baseline measures first.

### L3. Residual matching

**Meaning.** Match the leftovers of two partial decompositions instead of
requiring a full decomposition: states `s = (a, b, i_1..i_k)` with residual
`L(s) = aG + bQ − Σ P_{i_j}`; `L(s) = L(s')` is a relation over `F` (the
large-prime variation's "match the large primes"). **Status:** a toy module
(`residual_walk.rs`, strategies A, B, C1, C2 against plain rho R) exists;
it is not a framework oracle or collector and has never been measured on a
Koblitz base. **Toggle:** a new `relations.collector` (`residual`) in the IC
framework carrying the residual table; declared change is the collector only.
**Prediction:** work per relation `≈ √(r/B)`-type collision cost against the
direct oracle's `r/(poly(B))`-type scan cost, so at small `B` it wins by the
ratio of those two and loses once the direct hit rate exceeds the collision
rate; the exponent in `r` is unchanged (it is rho on the quotient). **Falsifier:**
no rung where relations/operation exceeds the baseline's. **Class:** at best
engineering; a win here does not move `B·D²/N`.

### L4. Orbit-invariant factor base

**Meaning.** A base `F` with `π(F) = F` and `−F = F`, so the fold of L1 is
exact. Two constructions: the Frobenius-invariant **subspace** (roots of a
linearised polynomial, needs `ord_n(2) < n − 1`), and the **orbit union**
(closure of a raw abscissa scan, `compact-orbit-scan`), the only one at
`n = 53, 61, 83, 131`. **Status:** both implemented. **Toggle:**
`koblitz-orbit:divisor=1` vs `binary-subspace` (prefix subspace of equal
dimension, not invariant, so orbits straddle the base and the fold is partial)
vs `compact-orbit-scan` at matched point count. **Prediction:** at equal
`|F|`, the invariant subspace and the orbit union have equal `columns`
(`|F|/2n`) and equal hit rate; the prefix base has `columns ≈ |F|/2` and a
lower fold. The subspace has an algebraic membership test (the lead for
algebraic oracles); the orbit union does not. **Falsifier:** hit rate per
column differing by more than the binomial interval between the two
invariant constructions (that would mean the subspace's points are not
uniform in the subgroup). **Class:** the fold is engineering; the algebraic
membership test is what L5/L7 need.

### L5. Phase-aware summation polynomials (ordered τ-slots)

**Meaning.** The "phase" of a summand is its Frobenius rotation `k` in
`τ^k(V)`. Decompose `R` with summand `i` constrained to `V^{(2^{j_i})}`, i.e.
solve `S_{m+1}(x_1^{2^{j_1}}, …, x_m^{2^{j_m}}, x_R)` for distinct phases
`(0, j_2, …, j_m)`. When `V` is **not** Frobenius-stable (every rung without a
subspace, hence `n = 131`), the `m` slots are different sets, the `S_m`
symmetry of the plain system is broken, and each solved system succeeds with
probability `|F|^m/N` instead of `|F|^m/(m!N)`; the relation still folds to
`V`'s columns with weights `λ^{j_i}`. **Status:** measured with an exact
oracle at `n = 11, 19, 23` (ratio 6.1 to 6.4 against `3! = 6`;
`ecdlp-hardness-work/ic-conductor-leads/LEADS.md` §3.1); **not implemented
in any production solver** (E1 result 3 in
`EXPERIMENTS_ISOGENY_CONDUCTOR_GAP_20261007.md`). On a Frobenius-stable `V`
the slots coincide and the lead is void, so it is measured on a prefix
(non-stable) base. **Toggle:** oracle parameter `slots=0;j2;j3` on the
algebraic and pair-table oracles. **Prediction:** trials per relation
`÷ m!` per solved system (`÷ 6` at `m = 3`), solver cost per system
unchanged, `n^{m−1}` systems available per target, LA unchanged.
**Falsifier:** a realised factor below 2 under the production solver (then
the solver was already exploiting the symmetry, as a symmetrised system
would). **Class:** engineering (`log₂ m!` bits of relation collection; it is
the *only* τ-dependent IC resource a non-stable base keeps, per LEADS V1).

### L6. Phase-aware quotient summation

**Meaning.** L5 applied to the **quotient** summation polynomial: the
`T₂`-symmetrised system in `u = 1/(x+1)`, `w = u² + u` (`koblitz_symmetrised`,
half the degree), with the slots as Frobenius conjugates of the `w`-subspace.
Two symmetries interact: symmetrisation already removes part of the `S_m`
symmetry on the *plain* side, so `m!` is an **upper bound** on what slots add
on top of it. **Status:** symmetrised systems implemented (`koblitz-symmetrised`
base, `symmetrised` oracle, `m ≤ 3`, `n ≤ 24`); slots not. **Toggle:** the same
`slots` parameter on the `symmetrised` oracle. **Prediction:** the realised
slot factor on the symmetrised system is between 1 and `m!`, and the
deliverable is that number. **Falsifier:** none needed; this row is a
measurement. **Class:** engineering.

### L7. Galois-orbit representations (equivariant and trace-reusing solving)

**Meaning.** Exploit the Frobenius `C_n` action on the barrel decomposition
system in the Gröbner step: block-diagonalise the Macaulay matrix by
characters, which at `n = 131` form one Galois orbit over `F_{2^130}` plus
the trivial character (G2); reuse the F4 trace across same-shape systems
(G1). **Status:** G2 scoped with one measured falsifier: the free `C_n`
action has the slice `k₀ = 0`, so gauge-fixing already recovers the `n`-fold
saving and the equivariant blocks then pay extension-field arithmetic
(`RESEARCH_G2_EQUIVARIANT_F4_SCOPE.md` §2); G1 unimplemented. **Toggle:**
solver engine (`matrix-f4` vs `matrix-f4:trace=replay`). **Prediction (G1):**
per-system F4 wall and word operations `÷ 3` to `÷ 100` at fixed shape, same
basis on every target. **Prediction (G2):** no gain over gauge-fixing unless
combined with the `S_m` symmetry; record the negative. **Falsifier (G1):**
a replayed trace giving a different basis on any target. **Class:**
engineering; moves the algebraic frontier by rungs of `n`, not to 131.

### L8. Conductor gap between curves (isogenies)

**Meaning.** The ECC2K-130 isogeny class has `cond(Z[π]) = 263 · p`,
`p ≈ 2^57.3`; the class is `{K_0}`, 262 curves at the 263-level (reachable,
kernel over `F_{q²}`), and `≈ 2^65` curves below the `p`-level that no known
method can write down or reach. **Status:** structure computed and verified
twice; V7, V8, V10 closed; E0, E1 (T23 only), E4, E5, E6 done; E8(1) analysed.
Measured so far: on floor curves the only IC resource lost is `τ` (orbit
folding `2n` in relations, realised τ-slot factor 1.0 in the production
pipeline because slots are unimplemented); per-call solver cost is
level-blind. **Done in this session:** E7 as an arithmetic model
(`research/ic_leads_20261010/e7_vertical_step_costs.py`, table in that
directory's README): the two components of the `n = 131` class are separated
by `≥ 2^53` relative to rho under every route modelled, the 263-level costs
`2^13` to reach, and on the toy ladder the only floor levels whose descent is
cheaper than rho (the honest reachability criterion) are `n = 61/1951`,
`n = 83/6473`, `n = 97/751943`, marginally `n = 41/409`; `n = 71, 73` have
none, which is why they model the challenge. The `n = 23` floor curves are
rebuilt by the CM route (`e2_floor_curves_polclass.sage`) as a check against
the Vélu construction. **Open:** E2 (floor curves at the landed rungs
`n = 41, 61, 83` under the compact-orbit pipeline), E3 (C37 four levels),
E8(2) (odd-`N` Kani prototype at C37 level 73), E9 (equal-precompute rho
across levels). **What this ladder adds:**
the floor curve is the exact control for "no Frobenius": a lead that claims a
gain *from* `τ` must vanish on the floor curve, and a lead that does not
(canonicalisation cost, residual matching, trace reuse) must survive there.
Every lead row therefore has an optional **floor-curve twin** at `n = 23`
(conductor 967, curves built in `research/isogeny_conductor_gap_20261007/`).
**Prediction:** L1, L4, L5, L6 vanish on the floor; L2, L3, L7 do not.
**Class:** control, not a lead; a hardness *difference* is reported only
under the V11 pre-registration.

## 5. Budget and order

| step | arm | rungs | cost on this host (load ≈ 100) | depends on |
|---|---|---|---|---|
| 1 | baseline sweep B1–B4, with L1 and L4 toggles (all plug-ins exist) | 17, 23, 31, 41 | minutes to an hour, operation counts only | release build of `ic` |
| 2 | L2 canonicalisation share from the counted oracle; prototype key (a) | 23, 31 | an afternoon | step 1 |
| 3 | L5 slots on the pair-table oracle over a prefix base (exact count of the realised factor) | 17, 23 | a day of Rust | step 1 |
| 4 | L5/L6 slots on the algebraic oracles | 17, 23 | two days | step 3 |
| 5 | L8 E7 cost table (arithmetic only, no Sage) and E2 at `n = 41` | 23, 41 | hours; `n = 41` floor curve needs `polclass` | none / Sage |
| 6 | L3 residual collector in the framework | 17, 23, 31 | two days | step 1 |
| 7 | L7 G1 trace replay | 31 | two days | none |

## 6. Non-claims

No row here is an online speedup, a total-work crossover, or a statement
about the challenge curve. The ladder's rungs are at most 41-bit subgroups;
the challenge is 129 bits and has no Frobenius-stable subspace, so the
subspace rows do not extrapolate to it and the orbit-union arm is the one
that does. Wall time from this session is not evidence.

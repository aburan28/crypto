# Lead ladder, round 2

**Arity note:** all `mitm` rows here are 2-sums (`m` defaults to 2); the n = 41 arms with `m = 3` are in a later section. (2026-10-10/11): L2, L3, L5 (table form), L8 floor twins

Same protocol as `RESULTS.md`; raw reports under `runs/20261011T00*/` with
commands, commit and binary hash. Binary rebuilt from this branch with the
residual collector, the canonicalisation price, the orbit-union base in
`ic bench`, and the floor-curve instance (`--model-b`). All runs recovered
and verified their planted logarithm. Host load 30–40 after a reboot; wall
time still not evidence.

## L2. Canonicalisation: priced, and negligible on the table pipeline

`ns_per_canon` is now measured in the binary calibration (least rotation of
the normal-basis coordinates, with shift):

| rung | ns per canonicalisation | ns per group addition | canonicalisation / addition |
|---:|---:|---:|---:|
| n = 13 (`K_0`) | 13.0 | 205 | 6.3 % |
| n = 17 (`K_1`) | 17.3 | 233 | 7.4 % |
| n = 23 (`K_0`) | 4.6 | 238 | 1.9 % |

A table probe is one addition, one canonicalisation and one lookup, so the
canonicalisation's share of a probe is 2–7 % here, and its share of the
whole decomposition phase (which also pays the walk step) is below 4 %.
**The ceiling for any cheaper orbit-invariant key is therefore ≈ 1.04× on
this pipeline**, and smaller on the compact-orbit pipeline, whose probe is
an S3 root solve. L2 is closed as measured: not worth a prototype key.
(The n = 23 figure is lower because the 23-bit coordinates fit three byte
tables; the cost is `O(⌈n/8⌉ + n)` word operations.)

## L3. Residual matching: implemented, measured, negative as the birthday bound predicts

`targets = residual` (`ic_framework::residual`): an r-adding walk that
subtracts one hash-chosen base point per step and stores every residual;
equal residuals (up to sign) give a relation over the base, a residual in
`±F` gives a classic relation. Same base and table as the baseline (the
table is built and unused, so `setup` is identical; read the relation
phase):

| rung | arm | trials | relations | trials / relation | decomposition adds | LA row ops | relation-phase ratio (adds) |
|---:|---|---:|---:|---:|---:|---:|---:|
| n = 13 | walk + `mitm-frobenius` | 16 | 16 | 1.0 | 200 | 208 | 1 |
| n = 13 | residual | 57.5 | 57.5 | 1.0 (all direct hits: the base is half the group) | 796 | 516 | 4.0× |
| n = 17 | walk + `mitm-frobenius` | 18.5 | 6 | 2.6 | 283 | 61 | 1 |
| n = 17 | residual | 808 | 8 | 101 | 1,960 | 258 | 6.9× |
| n = 23 | walk + `mitm-frobenius` | 62.5 | 27.5 | 2.6 | 430 | 268 | 1 |
| n = 23 | residual | 10,700 | 46 | 233 | 22,600 | 34,600 | 53× |

Readings. A residual collision is a birthday event on the subgroup: at
n = 23 the walk needed 233 steps per relation against `√(πr/2) ≈ 1,800` for
a cold rho solve, cheaper than one rho because the table keeps every
residual of every walk and later walks collide with it, but far above the
2.6 trials per relation of the direct 3-sum table. The gap widens with `n`
(4×, 7×, 53×) exactly as the comparison `√r` against `6r/|F|²` says it
must once `|F| > 2.4·r^{1/4}`. The relations are also denser (46 rows and
34,600 row operations against 27.5 and 268: each residual row carries a
walk's worth of subtracted points), so the LA loses too. **L3 is closed as a
negative on this family of bases**; the only regime where it could win is a
base below `r^{1/4}` points, which no rung uses. Classification:
engineering, and it does not move `B·D²/N` (1.2·10³ → 6.2·10⁵ at n = 23).

## L5. Ordered τ-slots, table form: the orbit-union base at equal points and columns (n = 23)

`compact-orbit-scan:columns=45` builds the signed-Frobenius **closure** of
the first 45 projected abscissae (2,070 points, 46 per column) — the base
shape that exists at `n = 131` — against the invariant subspace of the same
column count (2,025 points) and the prefix subspace of the same point count
(2,079 points, 1,040 columns). In the closure every rotation `τ^t P` of a
seed is a usable summand with weight `λ^t` on `P`'s column, which is
exactly the ordered-slot relation; the Frobenius-folded table probes all
rotations at once.

| base | columns | trials per relation | relations needed | table adds | base-build GAE (measured price) |
|---|---:|---:|---:|---:|---:|
| invariant subspace (`koblitz-orbit`) | 45 | 2.61 | 27.5 | 47,565 | 182 |
| orbit union (`compact-orbit-scan`) | 45 | **1.53** | 30.5 | 47,565 | 2,870 |
| prefix subspace (`binary-subspace`), negation only | 1,040 | 4.28 | 362 | 1,081,000 | 113 |

Readings. (1) The orbit union decomposes a target **1.7× more often per
trial** than the invariant subspace at equal points and columns, with the
same table cost; the relations needed are 11 % higher (orbit rows are
denser). The invariant subspace is not the better base; its only advantage
is the algebraic membership test the algebraic oracles need. (2) Against the
prefix base of equal point count the closure wins **2.8× in trials per
relation and 23× in columns**: that is the realised τ-slot gain in table
form, comparable to the `m! = 6` per-system figure of the exact-oracle
measurement once the fold (`n`) is separated from it (2.8 ≈ 6 × 0.47; the
prefix base's own fold-free yield is the denominator). (3) The base build
is 16× dearer in measured GAE because the scan's lifts and projections are
priced; it is target-independent and the one-time cost is 2,870 against a
relation phase of ≈ 1,100 per target.

So the τ-slot lead **is already realised in the table pipeline** by the
orbit-union base; what remains unimplemented is the algebraic-oracle form
(S3/S4 with per-summand rotations), which only matters when the pair
table does not fit, and that is the `n = 131` situation
(`RESEARCH_ECC2K130_IC_FEASIBILITY.md` §2).

## L8. Floor twin at n = 23: level-blind, as pre-registered

The same prefix base (dimension 11) and negation-folded table on `K_0` and
on the floor curve `b = 0x6a6f43` (orbit O00 of the 967-level, built by the
CM route, no Frobenius endomorphism; `ic bench --model-b`):

| curve | points | columns | trials per relation (8 repeats) | table adds | rho reference (`A = 2`) |
|---|---:|---:|---|---:|---:|
| `K_0` | 2,079 | 1,040 | 3.9–4.6, median 4.28 | 1,081,000 | S = 1.16 |
| floor O00 | 2,073 | 1,040 | 4.1–4.6, median 4.42 | 1,074,000 | S = 1.02 |

Within the spread of eight repeats. Everything the Koblitz model buys is in
the other rows: the `2n` fold (L1) and the closure's yield (L5), both of
which need `τ` and therefore vanish here, as predicted in the registry.
The hardness-difference pre-registration (V11) is unmet: no difference
beyond the Frobenius factors.

## Not reached in this round

- L5/L6 algebraic form (per-summand rotations in the S3/S4 systems and in
  the symmetrised system): not implemented.
- L7 G1 trace replay: not implemented.
- E2 at n = 61 (floor curves built, 1,952 of them) and n = 83: the
  compact-orbit pipeline (`koblitz_rank_fixture`) has no floor-curve input;
  `ic bench --model-b` covers the table pipeline only.

## Opportunity 1 (methodology §7): sizing the orbit-union base to the family optimum, n = 23

`compact-orbit-scan` at 5, 10, 20, 45 columns (230 to 2,070 signed points),
Frobenius-folded table, 8 targets each, all verified:

| columns | points | trials per relation | relations | table adds | decomposition adds | total GAE | `S` | `S` / rho (`A = 46`, S = 1.06) | `B·D²/N` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 230 | 59.1 | 5 | 690 | 647 | 2,370 | **1.64** | **1.55** | 204 |
| 10 | 460 | 16.6 | 9 | 2,530 | 505 | 4,360 | 3.01 | 2.85 | 312 |
| 20 | 920 | 4.84 | 16.5 | 9,660 | 438 | 12,000 | 8.32 | 7.9 | 549 |
| 45 | 2,070 | 1.53 | 30.5 | 47,565 | 420 | 51,600 | 35.6 | 33.8 | 1,190 |

Readings. The table cost falls with the square of the base, the relation
phase barely moves (647 against 420 additions: the walk step is most of a
trial and the probe loop is `|F|` additions per trial, so trials × `|F|`
is roughly constant), and the total drops **22×** from the 45-column base
to the 5-column base. At 5 columns the one-target total is within 1.55× of
the counted signed-Frobenius rho on the same instance; this is the family
optimum `|F| ≈ (#E·t/k)^{1/3} ≈ 200` of `ic_boundary::family_optimum_base`
landing where the formula says. It is a constant-factor result on a
21-bit subgroup and it is **not** sub-rho (`B·D²/N = 204`), but it means
every subspace-ladder row above was 10–20× oversized for a one-target
total, and the leaderboard's Koblitz rows should be re-based at the
optimum before any lead is read against them. The same sweep at n = 41
(columns 20–600) is in `runs/` once it completes.

## E2 inputs built this round

Floor curves by the CM route, all with `#E = #K_0` on the `a₂ = 0` twist:
n = 23 (966 curves, 42 orbits; identical to the Vélu family), n = 41
(410 curves, 10 orbits; `e2_floor_n41_ell409.json.gz`), n = 61 (1,952 curves,
32 orbits of 61; `e2_floor_n61_ell1951.json.gz`). The n = 41 twin under the
table pipeline is in `runs/20261011T001107Z`.

## Base sizing at n = 41 (K_0, 39-bit r, m = 2, 4 targets each)

`compact-orbit-scan` at 20, 60, 200 columns (82 points per column), Frobenius-folded table:

| columns | points | trials per relation | relations | table adds | decomposition adds | total GAE | `S` | rho `S` (`A = 82`) | IC / rho |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 | 1,640 | 363,000 | 14.5 | 17,200 | 5.80·10⁶ | 6.32·10⁶ | 8.5 | 0.21 | 41× |
| 60 | 4,920 | 37,600 | 37 | 150,000 | 1.31·10⁶ | 1.58·10⁶ | **2.13** | 0.21 | 10× |
| 200 | 16,400 | 4,030 | 127 | 1.65·10⁶ | 5.01·10⁵ | 2.22·10⁶ | 2.99 | 0.21 | 14× |

The optimum for 2-sums at 39 bits sits near 60 columns, where the table and
the relation phase are within 10× of each other; the 600-column ledger
shape is on the table-dominated side. Against the matched rho the best
2-sum row is 10× behind, with the gap opening from 1.55× at n = 23: the
`r^{1/6}` law of the pair-table family (`family_optimum_s`), not a lead.

## n = 23 with three summands (the ledger's arity)

Every target decomposes (trials per relation 1.0 on all three bases at
2,000 points, 21 bits), so the arity moves the comparison from "how often"
to "how many probes": the 3-sum oracle probes `R − P` for base points
until a hit (52–91 lookups per relation on the folded bases, 2,790 on the
prefix base). The ratios of the 2-sum rows carry over: orbit union against
invariant subspace 0.89 relations, 0.58 LA row operations, 0.92
decomposition additions; prefix base 18× the relations and 8.8× the
decomposition additions of the baseline. Frobenius-folded against
negation-folded table: 21.6× on setup, nothing else.

## Base sizing at n = 41 with three summands (4 targets each)

| columns | points | trials per relation | relations | probes (lookups) | decomposition adds | total GAE | `S` | IC / rho (`A = 82`, S = 0.21) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 20 | 1,640 | 898 | 20 | 2.86·10⁷ | 2.86·10⁷ | 3.10·10⁷ | 41.7 | 199× |
| 60 | 4,920 | 27.3 | 55 | 7.24·10⁶ | 7.24·10⁶ | 7.90·10⁶ | 10.8 | 51× |
| 200 | 16,400 | 1.35 | 186 | 1.59·10⁶ | 1.59·10⁶ | 3.36·10⁶ | 4.57 | 22× |

Against the 2-sum rows at the same columns (8.5, 2.13, 2.99): **three
summands lose to two at every size tried**, because the 3-sum oracle pays
`|F|` probes (one addition each) per trial and the extra yield does not
repay them until the base is large enough that trials per relation reach
1; at 200 columns it is 1.35 and still 2× behind the 2-sum optimum. The
600-column `m = 3` run (the ledger's shape) is in `runs/20261011T001608Z`
when it completes. For the pair-table family at 39 bits the operation-
counted optimum is therefore the 2-sum table near 60 columns; the
leaderboard's `m = 3` collection rows are priced in a different unit and
window and are not contradicted by this, but a re-based row should start
from the 2-sum optimum.

## L5, algebraic form: ordered τ-slots in the descended S3 (new this round)

New plug-ins: `koblitz-closure:dimension=d` (the signed-Frobenius closure
of the prefix subspace `V` of dimension `d`, one column per orbit, `V`'s
basis kept) and `descent-algebraic:slots=all` (per target, one descended
system per rotation `j` with summand 2 over `σ^j(V)`, first lifted hit
wins; `pq_descent_symbolic::descend_with_bases`). Solver `fes-f2`
(exhaustive over the `2d` Boolean unknowns, exact), so the comparison is
about yield, not solver heuristics. 8 targets each, all verified.

| rung | base | arm | targets tried | relations | systems solved | **hits per system** | relation-phase GAE |
|---:|---|---|---:|---:|---:|---:|---:|
| n = 17 (`K_1`) | closure of dim 6: 986 pts, 29 cols | plain descent (`V × V`) | 6,799 | 126 | 6,744 | 0.0187 | 4.56·10⁵ (median run) |
| n = 17 | same | slots = all (17 systems/target) | 382 | 133 | 5,313 | 0.0250 (**1.34×**) | 1.90·10⁵ |
| n = 23 (`K_0`) | closure of dim 7: 2,850 pts, 63 cols | plain descent | 174,614 | 210 | 173,639 | 0.00121 | 9.14·10⁶ |
| n = 23 | same | slots = all (23 systems/target) | 5,972 | 290 | 134,174 | 0.00216 (**1.79×**) | 4.69·10⁶ |
| n = 23 | same | `mitm-frobenius` table (reference) | 564 | 232 | — | — | 1.09·10³ + 9.3·10⁴ setup |

Readings. (1) The pre-registered prediction was `m! = 2` hits per solved
system; measured **1.79 at n = 23** (binomial 95 % interval on 290 hits
≈ ±0.21, consistent with 2) and **1.34 at n = 17** (`K_1`, cofactor 2,
986 points in a 2^16 group: the coset structure LEADS saw inflate the
n = 11 figure deflates this one). (2) Because every rotation is a separate
system, the slots also multiply the relations **per target** by ≈ `2n`
(trials per relation 18.3 against 897 at n = 23), which is the relation
space the table pipeline gets for free; the solver cost per relation halves
(4.7·10⁶ against 9.1·10⁶ GAE). (3) At these sizes the algebraic oracle is
50–100× dearer than the pair table on the same base; the lead matters where
the table does not fit. The plain arm hits only at rotation 0 by
construction; the slotted arm's hits are spread over all rotations
(`slot_hit_rotation_j` counters in the reports), which is the structural
check that the per-summand bases are doing what they claim.

## A1. `Z[τ]`-multiplier closure (new `norm` parameter of `compact-orbit-scan`), n = 23, first run

Multipliers `α = a + bτ` of norm `a² − ab + 2b² ≤ norm`; `α(P) = [a]P + [b]τ(P)`
joins `P`'s column with weight `a + bλ`, its whole signed Frobenius orbit
with it. 5 columns, Frobenius-folded 2-sum table, 8 targets, all verified:

| norm | points | points per column | trials per relation | table adds | `S` |
|---:|---:|---:|---:|---:|---:|
| 0 (plain `±τ^k`) | 230 | 46 | 59.1 | 690 | 1.63 |
| 4 | 690 | 138 | 27.6 | 2,070 | 3.60 |
| 8 | 1,150 | 230 | 18.3 | 3,450 | 5.69 |
| 16 | 2,300 | 460 | 9.4 | 6,900 | 10.4 |
| 45 columns, norm 0 (equal points to norm 16) | 2,070 | 46 | 1.53 | 47,565 | 35.6 |

Reading of the first run: the relations are correct (the `a + bλ` weights
verify), but the yield falls **linearly** with the points (3×, 5×, 10×
points → 2.1×, 3.2×, 6.3× fewer trials) where a 2-sum table over
independent points falls with the square, and the table grew linearly too
(690 → 2,070 → 3,450 → 6,900 additions). The cause is in the table, not the
lead: `FrobeniusPairTable::build` took **one representative per column**
and paired it against the base, which is complete only when a column is a
single orbit; with `M` orbits per column it builds `M·K × |B|` pairs instead
of `|B|²/(4n)`. Fixed (one representative per *orbit*, found by the
canonical abscissa); the rerun follows below.

## L8 floor twin at n = 41, three summands (4 targets)

Floor curve `b = 0xbe5ad65b0b` (orbit O00 of the 409-level), prefix
subspace of dimension 13 (8,060 points, 4,030 columns), negation-folded
3-sum table: trials per relation **26**, 3,700 relations, 7.6·10⁸ probe
additions, `S = 1,040` against the floor's negation-only rho `S = 1.02`
(the floor loses the `√(2n) = 9.1` of the signed-Frobenius walk as well).
The matching `K_0` row with the same base and oracle is the last arm of
`runs/20261011T001608Z` and is compared here when it lands; the
level-blindness prediction is trials per relation within the spread of 26.

### A1 rerun with the per-orbit pair table (n = 23, 5 columns, 2-sums, 8 targets, all verified)

| norm | points | points per column | trials per relation | relations | LA row ops | table adds | base build GAE | total GAE | `S` |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | 230 | 46 | 59.1 | 5 | 46 | 690 | 319 | 2,360 | 1.63 |
| 4 | 690 | 138 | 10.9 | 5 | 45 | 5,520 | 1,980 | 8,590 | 5.93 |
| 8 | 1,150 | 230 | 4.83 | 5 | 41 | 15,000 | 3,670 | 19,700 | 13.6 |
| 16 | 2,300 | 460 | **1.50** | 6 | 47.5 | 58,600 | 7,080 | 66,800 | 46.1 |
| 45 columns, norm 0 | 2,070 | 46 | 1.53 | 30.5 | 382 | 47,565 | 2,860 | 51,600 | 35.6 |

Readings. (1) With the table complete, the yield scales with the **square**
of the points as the lead claims (3×, 5×, 10× points → 5.4×, 12×, 39× fewer
trials, saturating near one trial per relation), and the table grows with
the square too (58,600 at 2,300 points against 47,565 at 2,070): the
multiplier closure behaves as a set of independent points for the 2-sum
table. (2) At equal points (norm 16 against 45 plain columns) the yield is
equal (1.50 vs 1.53), the table is equal, and the columns are **9× fewer**
(5 vs 45): relations needed 6 against 30.5, LA 8× smaller. (3) The price is
the base build, 7,080 against 2,860 GAE, because every `α(P)` costs two
scalar multiplications, which at this size outweighs the LA saving, so the
one-target total is 1.3× worse. (4) The lead therefore pays exactly where
the feasibility note says the compact-orbit pipeline is bound: when
relations are expensive (probes per relation ≫ 1) and the rank stage needs
`K + 1` of them, `K ÷ M` at equal table and equal yield is a direct `M×` on
the rank precompute. On the toy table at 21 bits the relation phase is
already one trial per relation, so there is nothing left to save there.
Prediction for the compact-orbit pipeline at the landed rungs (K = 600,
probes per relation 10⁶–10⁸): rank work ÷ `M` at equal table memory, with
`M ≈ 10` at norm 16, i.e. the same `K` with a `100×` smaller root table or
the same table with `10×` fewer relations; unmeasured there.

## §7 items 2, 3, 6, 7: status

- **Re-basing the leaderboard (item 2).** The inputs are now measured: the
  orbit-union base with a 2-sum Frobenius-folded table at the family
  optimum (3–5 columns at n = 23, ≈ 60 at n = 41) is the operation-counted
  optimum of the pair-table family on `K_0`, 1.5× behind the counted
  signed-Frobenius rho at 21 bits and 10× behind at 39 bits; three summands
  lose at both rungs. Table A of `LEADERBOARD.md` is priced in a different
  unit (batched additions per `√r`, `ic price`) and window (one target
  after collection), so the re-based rows have to be produced by
  `ic price` on `compact-orbit-scan` bases of the optimum size, which is a
  measurement campaign, not a code change; the ladder here fixes the base
  shape and size that campaign should use.
- **Charge-free base build (item 3).** After the A1 measurement the build
  cost that matters is the multiplier images (`7,080` GAE at norm 16
  against `2,860` plain), not the inversions of the plain scan; the
  multipliers are 2-bit scalars whose double-and-add is already 2–3 group
  operations, so an addition chain over the shared `P`, `τ(P)` would save
  at most a third of it. Not pursued; the plain scan's `_uncharged` lifts
  remain as recorded.
- **Batched probes (item 6).** Invisible in operation counts and only
  measurable in the wall-clock window on an isolated host; not pursued on
  this host.
- **Calibration pins (item 7).** the `icv1-f2m23-t5197-69e76b73` entry added to
  `docs/ic/calibration.json` from the median of 10 reports (provenance in
  `sources`); `K_0` n = 41 was already pinned, which is why its reports
  show `pinned` units. The pins are from this loaded host and should be
  re-measured on the ledger's isolated host before any `S` from them is
  compared across hosts.

## Next-step scoping (2026-10-11): L6 and L7, decided by reading the code

- **L6, symmetrised slots.** `build_symmetrised_system` takes one `V ∋ 1`
  for all summands (`w_basis = AS(v_basis[1..])`); the slot form needs one
  `(v, w)` pair per summand with `w_i = AS(σ^{j_i} v)`, which commutes with
  `σ`, a `u`-frame closure base that keeps the `T`-translation, and a lift
  with per-summand bases: ≈ 150 lines mirroring `descend_with_bases`. The
  expected answer is what decides whether to spend them: the symmetrised
  variables `(w, s)` are already the symmetric functions of the summands,
  so the plain symmetrised system counts an unordered decomposition once;
  the `m!` the slots recover on the plain `x`-system is therefore mostly
  gone before the slots start, and the realised factor should sit near 1
  (the registry's "between 1 and `m!`"). Deferred as low value; the table
  form of L5 already runs on `F_u`.
- **L7, F4 trace replay (G1).** `pq_f4_f2::groebner_basis_f4` selects every
  pair of minimal degree per step (normal strategy) and keeps per-step
  statistics but not the identity of the rows that reduced to zero. The
  Joux–Vitse record is, per step, the `(multiplier, polynomial)` keys of the
  rows that survived; a replay mode would form only those rows and skip
  symbolic preprocessing of the rest. Polynomial indices are target-
  dependent in general, so the trace must be keyed by position in the
  degree-sorted basis, which the same-shape argument makes stable. Half a
  day of Rust with a correctness gate (same reduced basis on every target);
  the measurement is the n = 31 F4 cell. Not started.

### A1 where relations are dear: n = 41, three summands, equal points (4 targets each, all verified)

| base | points | columns | trials per relation | relations | probes (adds) | table adds | base build GAE | total GAE | `S` | IC / rho (`A = 82`) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 60 plain columns | 4,920 | 60 | 27.3 | 55 | 7.24·10⁶ | 1.5·10⁵ | 6,830 | 7.90·10⁶ | 10.8 | 51× |
| 20 columns, norm 4 | 4,920 | 20 | 35.2 | 20.5 | 3.24·10⁶ | 1.5·10⁵ | 11,400 | 3.50·10⁶ | 4.94 | 23.6× |
| 6 columns, norm 16 | 4,920 | 6 | 21.5 | 6.5 | 6.55·10⁵ | 1.5·10⁵ | 11,700 | 8.71·10⁵ | **1.17** | **5.6×** |

The prediction held: equal points give equal yield (21.5–35.2 trials per
relation, within the spread of 4 targets) and an identical table (the
Frobenius-folded table is a function of the point set), while the
relations needed fall with the columns (55 → 6.5), so the probe count,
which is the dominant phase at `m = 3`, falls 11× and the one-target total
**9.2×**. The base build is 1.7× dearer (two scalar multiplications per
image) and is now 1.3 % of the total. This is the first row on the ladder
within one decimal order of the matched rho at 39 bits, and it is the
constant-factor mechanism the feasibility note's rank-stage accounting
(`r/(n²K)` probes with `K` columns) says the compact-orbit pipeline needs:
`K ÷ 10` at equal table memory, i.e. a 10× smaller rank precompute, or the
same `K` with the yield of a 10× larger table. Still `B·D²/N ≫ 1` and still
a pair-table constant, not an exponent.

### A1 limit at n = 41, three summands: fewer columns, larger norm (4 targets each, all verified)

| base | points | columns | trials per relation | relations | probes (adds) | table adds | base build GAE | total GAE | `S` | IC / rho |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 6 columns, norm 16 | 4,920 | 6 | 21.5 | 6.5 | 6.55·10⁵ | 1.5·10⁵ | 11,700 | 8.71·10⁵ | 1.17 | 5.6× |
| 3 columns, norm 32 | 4,920 | 3 | 34.8 | 4 | 6.69·10⁵ | 1.5·10⁵ | 11,600 | 8.86·10⁵ | 1.19 | 5.7× |
| 2 columns, norm 64 | 6,230 | 2 | 22.5 | 3 | 4.05·10⁵ | 2.4·10⁵ | 15,200 | 6.94·10⁵ | **0.94** | **4.5×** |
| 1 column, norm 160 | 7,790 | 1 | 37.8 | 2 | 5.75·10⁵ | 3.74·10⁵ | 18,900 | 1.02·10⁶ | 1.37 | 6.5× |

The floor is in sight: a one-column base still needs two relations (the
target column and one more), each costing `trials × |B|` probes, so below
three columns the relation phase stops falling and the table (∝ `|B|²/4n`)
takes over. The operation-counted optimum of the pair-table family on `K_0`
at 39 bits with everything this ladder has (orbit union, `Z[τ]` closure,
Frobenius-folded table, three summands) is **≈ 2 columns of ≈ 6,000 points,
`S ≈ 0.94`, 4.5× the matched signed-Frobenius rho**, against 51× for the
ledger-shaped 60-column base and 199× for 20 plain columns. Still a constant
on a toy rung; the exponent is the family's.

### L5 algebraic form, three summands (descended S4, `exhaustive` solver on 12 unknowns)

Only summand 2 is rotated (`slots=all`), summands 1 and 3 stay over `V`, so
the plain system's `3! = 6` symmetry is broken down to the swap of summands
1 and 3: the prediction is **3** hits per solved system relative to plain.

| rung | base | arm | targets tried | relations | systems solved | hits per system | ratio |
|---:|---|---|---:|---:|---:|---:|---:|
| n = 17 (`K_1`) | closure of dim 4: 136 points, 4 columns | plain descended S4 | 22,982 | 35 | 22,451 | 0.00156 | 1 |
| n = 17 | same | slots = all (17 systems per target) | 686 | 34 | 11,387 | 0.00299 | **1.92** |
| n = 23 (`K_0`) | closure of dim 4: 370 points, 9 columns | plain descended S4 | 394,971 | 69 | 389,587 | 0.000177 | 1 |
| n = 23 | same | slots = all (23 systems per target) | 8,197 | 71 | 187,698 | 0.000378 | **2.14** |

The measured factor is about two thirds of the prediction at both rungs
(1.92 at n = 17, 2.14 at n = 23, against 3), the same shortfall as at
`m = 2` (1.34 and 1.79 against 2). With 69–71 hits per arm the binomial
interval on the n = 23 ratio is roughly ±0.5, so 3 is not excluded there;
the consistent shortfall across four rows suggests a systematic cause
(the closure's dimension-4 seed subspace overlaps its own Frobenius
images, so some rotated systems duplicate plain ones) rather than noise.
Relations per target rise ≈ 48× at n = 23 (trials per relation 115 against
5,724) because every rotation is its own system. Lift failures (18 and 10) are abscissae whose points lie on
the twist, as the oracle documents; every recovered logarithm verified.

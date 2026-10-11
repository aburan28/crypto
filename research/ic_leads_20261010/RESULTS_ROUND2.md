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
(410 curves, 10 orbits; `e2_floor_n41_ell409.json`), n = 61 (1,952 curves,
32 orbits of 61; `e2_floor_n61_ell1951.json`). The n = 41 twin under the
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

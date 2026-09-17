# A different decomposition solver, run on the real curve

**Experiments:** `research/nagao_relations/solver_10`, `research/nagao_relations/solver_11`
**Frozen contracts:** `solver_10/contract.json`, `solver_11/contract.json` (both committed before execution)
**Evidence:** `solver_10/raw.jsonl` (128 trials), `solver_11/raw.jsonl`
**Background:** [`RESEARCH_ECC2K130_DECOMPOSITION.md`](RESEARCH_ECC2K130_DECOMPOSITION.md)
(existence, admissibility, and the `m·2^131` oracle bound),
[`research/nagao_relations/README.md`](research/nagao_relations/README.md)
(the Riemann–Roch encoding),
[`research/nagao_relations/scaling_23_29.md`](research/nagao_relations/scaling_23_29.md)
(the same solver family on the toy ladder)

**The question.**  The Riemann–Roch (Nagao) decomposition solvers in
`solver_05/07/08` beat the Semaev controls on the toy ladder, `n ≤ 29`.  Toy
wins are cheap.  Does the advantage survive on `K_0 : y² + xy = x³ + 1` over
`F_2^131` — the ECC2K-130 curve itself, group order checked against the Koblitz
recurrence as `4r` with the published 129-bit prime — with the factor base,
arity, targets, modes and budget held identical across every variant?

**Bottom line.  Yes, it survives, and it is not close.  It also does not
matter, and §4 says exactly why in the attack's own units.**  Those are two
separate results and conflating them is how this route gets over-sold:

- **The solver comparison (§2).**  At `d = 6`, the two RR solvers resolved
  **32 of 32** matched slots; both Semaev controls resolved **0 of 32** under
  the same 120 s all-phase budget.  Across the whole campaign, **128 trials,
  zero errors, 28 independently verified three-summand relations on the real
  curve**, every complete planted run recovering its planted triple.
- **What it is worth (§4).**  Nothing about the attack.  The panel prices an
  *oracle* at `d ∈ {6, 7}` against an existence threshold of `d = 45`.  A
  uniform target at `d = 6` decomposes with probability at most `5.4 × 10^-36`,
  and closing that gap needs `B ≈ 1.6 × 10^13` stored points.
- **The batch result, and its missing control (§5).**  Amortising setup across
  a batch of targets cannot be the source of an index-calculus advantage,
  because rho amortises better.  The RR solver's target-independent share is
  measured at **0.12 %** (`d = 6`) and **0.07 %** (`d = 7`).  Batching 16
  targets buys the RR oracle 0.1 %; it buys rho 75 %.
- **Two structural limits (§6),** read off the construction rather than
  measured: the arity is fixed at three, and the factor base must be an
  `F_2`-subspace.  "More summands" and "a Hamming-weight base" are outside this
  encoding, not merely untested by it.

## 0. What is matched, and what correctness means with no oracle

Every variant is imported **unmodified** from its own frozen folder —
`quadratic-optimized` from `solver_07`, `quadratic-image` from `solver_08`,
`s4-symmetric` from `solver_06`, `chained-s3` from the `solver_04` frontend.
Only the driver is new.  Shared across all four: the factor base
`F = { P : x(P) ∈ V \ {0} }` with `V` the span of the first `d` ONB basis
vectors, arity `m = 3`, the target list, both modes, and a 120 s all-phase
budget per instance.

The toy campaigns could lean on an exhaustive pair oracle over the whole group.
At `n = 131` no such oracle exists.  Rather than quietly weaken the reused SAT
drivers, the two assertions that need that oracle are stood down **explicitly**
— a permissive `expected` object that removes checks and never adds one — and
replaced by independent re-verification in the driver: every reported abscissa
triple is re-lifted with `pointFromX`, re-checked for subspace membership and
the exclusions, and its signed group sum recomputed and compared with the
target.  A `complete` run on a planted target that misses its planted triple is
a correctness failure.

Two strata, and the difference matters:

| stratum | how built | what it tests |
|---|---|---|
| `uniform` | `pointFromX` of a uniform random abscissa | the cost of *exhausting* a candidate space that contains nothing |
| `known_decomposable` | `R = P₁+P₂+P₃` from a random factor-base triple, triple recorded | the cost of *finding* a relation that provably exists |

The planted stratum is a constructed existence witness.  **It is not a natural
yield estimate and nothing below treats it as one.**

## 1. The boundary, stated before anything ran

```text
#E(F_2^131) = 4·r,  r = 680564733841876926932320129493409985129   (prime, 2^129.0000)
rho reference:      2^60.8090      S = 0.077430
existence threshold: m·l ≥ 131 + log₂ m!   →   l = 45 at m = 3
```

The panel runs at `d = 6` and `d = 7`.  That is `18` and `21` against `131`, so
**the uniform stratum is expected to contain no decomposition at all**, and its
runs measure exhaustion cost, not yield.  This was written into the contract
before execution, not discovered afterwards.

## 2. The matched panel

Resolved slots within the shared 120 s all-phase budget, out of four attempted;
`rel` counts independently verified relations.  Factor base: 27 points at
`d = 6`, 66 at `d = 7`.

### d = 6

| variant | uniform first | uniform enum | planted first | planted enum | rel | median resolved |
|---|:--:|:--:|:--:|:--:|--:|--:|
| **quadratic-image** | 4/4 | 4/4 | 4/4 | 4/4 | 8 | 20.5–32.2 s |
| **quadratic-optimized** | 4/4 | 4/4 | 4/4 | 4/4 | 8 | 27.7–43.4 s |
| s4-symmetric | 0/4 | 0/4 | 0/4 | 0/4 | 0 | — |
| chained-s3 | 0/4 | 0/4 | 0/4 | 0/4 | 0 | — |

**32 of 32 against 0 of 32.**  On planted targets the image solver reached its
first verified relation in **16.0–22.8 s**, the optimized solver in
21.2–30.9 s; neither Semaev control resolved a single slot in either stratum or
either mode.  On uniform targets both RR solvers *completed* — that is, they
exhaustively certified that no decomposition of the target exists over this
factor base — in 32 s and 43 s.

Counted field API operations on the matched complete `d = 6` enumerations, same
13,588–13,608 candidate functions on each side:

| variant | uniform ops | planted ops | ratio |
|---|--:|--:|--:|
| quadratic-optimized | 7,844,652 | 7,840,178 | 1.000 |
| quadratic-image | 6,483,976 | 6,484,110 | **0.827** |

The image-space support test carries its toy-ladder advantage to `n = 131`:
same candidate sequence, 17 % fewer field operations.

### d = 7

The candidate space grows about fourfold from `d = 6` to `d = 7`, and the
budget boundary lands inside the panel:

| variant | uniform (any mode) | planted first | planted enum | rel |
|---|:--:|:--:|:--:|--:|
| quadratic-image | 0/8 | **4/4** | 0/4 | 8 |
| quadratic-optimized | 0/8 | 2/4 | 0/4 | 4 |
| s4-symmetric | 0/8 | 0/4 | 0/4 | 0 |
| chained-s3 | 0/8 | 0/4 | 0/4 | 0 |

The RR solvers still find and verify relations where one exists — the image
solver on all four planted targets, the optimized solver on two — while no
control resolves anything.  The `enum` rows are `0/4` **with relations found**:
the solvers verified relations and then ran out of budget before exhausting the
space, so the run is recorded as `timeout` and its solution set is a lower
bound, not a complete set.

**Timeouts are recorded as unknown, never as negative mathematical evidence.**
A `timeout` here bounds this implementation, this budget and this panel.  It
says nothing about whether a decomposition exists, and §6 of `AGENTS.md`
forbids reading it as if it did.

## 3. Correctness

128 trials, **zero errors**, **28 independently verified relations** on the
ECC2K-130 curve.  Every one was re-derived from its abscissas alone — lifted,
re-checked against the factor base and the exclusions, and re-summed in the
group — without trusting the solver that produced it.  All eight complete runs
on planted targets contained their planted triple.

## 4. What this is worth, in the attack's own units

**Nothing, and the arithmetic is not close.**  With a stored base of `B` points
and three summands, at most `C(B+2, 3)` distinct points are representable, so a
uniformly sampled target of the order-`r` subgroup decomposes with probability
at most `min(1, C(B+2,3)/r)`:

| `d` | `B` | `P(decompose) ≤` |
|--:|--:|--:|
| 6 | 27 | `5.37 × 10^-36` |
| 7 | 66 | `7.36 × 10^-35` |

Full coverage needs `B ≈ 1.6 × 10^13 ≈ 2^43.9` stored points — consistent with
the `2^44.5` in §6 of the background note.  The panel is therefore a
**solver-quality measurement at a base size the attack cannot use**, and the
background note's `m·2^131` bound is untouched by it: that bound is a statement
about the product `relations × targets × oracle`, and a faster oracle at fixed
`(m, l)` moves one factor while the law holds the product constant.

No `Λ`, no full-DLP `S`, and no rho ratio is computed from these runs.  Not
because they would be unflattering, but because at `d ≪ 45` they would be
meaningless.

## 5. Batching, and the baseline that decides it

A parallel line of work reported that sharing index-calculus setup across 16
targets used about 28 % fewer instructions and 22 % less time than **running
the rho solver separately for each target**, while a single target cost 7.8 %
more than rho.  That baseline is the thing under test.

Rho amortises across a batch too.  By the law this repository already enforces
in `src/ecc_safety.rs` (`check_multi_target_margin`, Galbraith–Lin–Scott /
Kuhn–Struik), `k` discrete logarithms in one group with one shared generator
cost about `√(k·r)` in total, so **per-target rho cost falls as `1/√k`, without
a floor.**

`solver_11` measures how much of the RR solver's work can amortise at all.  The
image spaces of `T_u(w) = w² + uw` on `V` depend only on `V`, so they are
genuinely target-independent; they were re-implemented independently here and
asserted equal to `solver_08`'s frozen construction before any cost was split.
Everything else — the geometry, the batch inversions, the search — depends on
the target.

| `d` | target-independent ops | per-target ops | **amortisable share** |
|--:|--:|--:|--:|
| 6 | 1,982 | 1,618,400 | **0.12 %** |
| 7 | 5,123 | 6,931,952 | **0.07 %** |

Per-target cost as a fraction of its `k = 1` value:

| `k` | RR oracle | amortised rho |
|--:|--:|--:|
| 1 | 1.000 | 1.000 |
| 16 | 0.999 | **0.250** |
| 1024 | 0.999 | **0.031** |
| 2^20 | 0.999 | **0.001** |

**Batching 16 targets buys the RR oracle one part in a thousand; it buys rho a
factor of four.**  The RR per-target cost has a strictly positive floor and is
already within 0.1 % of it at `k = 2`.  So for every `k`, the rho side gains at
least as much from the batch as the index-calculus side does: a batch of size
`k` cannot create a crossover that does not already exist at `k = 1`.

A batch win over `k` independent rho runs measures the baseline, not the
algorithm.  Whatever produced the 28 % figure, it cannot have been amortised
decomposition setup at this scale, and the correct comparison — batched against
batched — moves *against* index calculus as `k` grows.

## 6. Two structural limits

Read off the construction, not measured, and therefore not a matter of budget:

1. **The arity is fixed at three.**  The encoding lives in `L(4O)` with
   `f = x² + ax + c + by`; `div(f) = (P₁)+(P₂)+(P₃)+(−R)−4(O)` exhausts pole
   order four.  Adding summands is not a parameter of this solver — it requires
   a different Riemann–Roch space.
2. **The factor base must be an `F_2`-subspace.**  Both support tests are
   subspace tests: `solver_05` divides by the subspace polynomial `L_V`,
   `solver_08` tests membership in the image of the `F_2`-linear map
   `w ↦ w² + uw` on `V`.  A normal-basis Hamming-weight factor base is not a
   subspace, so it is outside this encoding entirely.

Both are worth stating because both are natural next knobs to reach for, and
neither is reachable from here.  Widening the base along a *different* axis —
higher `d`, which stays inside the encoding — is what §2's `d = 7` rows price,
and the cost grows about fourfold per dimension.

## 7. What this does and does not license

Established:

- The RR solver family's advantage over Semaev is **not an artefact of toy
  field sizes**; it reproduces at `n = 131` with a larger margin than on the
  ladder (32/32 against 0/32 at `d = 6`).
- `quadratic-image` is the strongest member at full size, on both resolved
  slots and counted operations (0.827× `quadratic-optimized`).
- Three-summand decompositions of real ECC2K-130 points are findable in
  seconds *when the base is small enough to search and the target is known to
  decompose* — 28 of them, each independently verified.

Not established, and not claimed:

- Any relation yield, any `Λ`, any `S`, any ratio to rho, any attack cost.
- Anything about the natural (uniform) case: at `d ∈ {6, 7}` it is
  overwhelmingly likely there is nothing to find, and the runs confirm the cost
  of establishing that, not a failure to find.
- Any crossover from batching, which §5 argues is the wrong place to look.

The open question this leaves is unchanged and is the one the background note
already names: an oracle that beats exhaustive search over its own candidate
set by `2^{70.19 + log₂ m}`.  Being 17 % cheaper per candidate, and beating
Semaev on every matched slot, is not that — it is a better constant on a
candidate count that the product law holds fixed.

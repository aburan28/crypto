# Target: make the Koblitz decomposition oracle reach a useful `m`

**Harness:** `src/cryptanalysis/koblitz_bench.rs`
**Baseline run:** `cargo run --release --example koblitz_scaling_bench`
**Machine-readable:** `… --example koblitz_scaling_bench -- --json`
**Under test:** `koblitz_index_calculus`, `koblitz_groebner`, `semaev_sat`
**Background:** `RESEARCH_KOBLITZ_INDEX_CALCULUS.md`

This is a scoped, measurable target with a single primary number, a
correctness gate that cannot be gamed, and five pre-registered
hypotheses.  It is written to be picked up as an autolab thread.

---

## The one-paragraph version

Index calculus on a Koblitz curve needs to write a target as a sum of
`m` factor-base points.  With a Frobenius-invariant factor base of size
`|F| ≈ 2^ℓ` and a subgroup of order `≈ 2^n`, the decomposition only
*exists* when `|F|^m / m! ≳ 2^n`, i.e. when `m ≈ n/ℓ`.  Every oracle we
have works at `m = 2` and `m = 3`.  None can reach `m ≈ n/ℓ` for any `n`
worth the name, and the reason is not the algebra — it is that chaining
`S₃` over `m − 2` intermediate points costs `(m − 2)·n` Boolean unknowns
on top of the `m·ℓ` that the factor base costs.  **Drive that number
down.**

---

## Primary metric

```
    unknowns(n, ℓ, m) = m·ℓ + (m − 2)·n          (koblitz_bench::n_vars_for)
```

evaluated at the **useful** summand count `m = ⌈n/ℓ⌉`, on the ladder of
`n` that have a non-trivial invariant subspace.  Lower is better.  The
current representation caps a solvable instance at 64 unknowns.

| n | ℓ | m = 2 | m = 3 | m = 4 | **m = ⌈n/ℓ⌉** |
|--:|--:|------:|------:|------:|--------------:|
| 5 | 4 | 8 | 17 | 26 | **8** |
| 7 | 3 | 6 | 16 | 26 | **16** |
| 9 | 6 | 12 | 27 | 42 | **12** |
| 15 | 4 | 8 | 27 | 46 | **46** |
| 21 | 6 | 12 | 39 | 66 | **66** |
| 31 | 5 | 10 | 46 | 82 | **190** |
| 63 | 6 | 12 | 81 | 150 | **633** |

Everything at or below 64 is reachable today; everything above is not.
The gap at `n = 31` is 3×, at `n = 63` it is 10×.

**Secondary metrics** (report all three, every run):

1. Largest `(n, m)` whose decomposition an oracle *solves* within
   budget, with all oracles agreeing.
2. Median ms per target, per oracle (`format_oracle_table`).
3. First fall degree of the system (`SystemProfile::fall_degree`).

---

## Correctness gate — non-negotiable

`InstanceBench::disagreements` **must be 0** in every reported run. It
counts two things: targets where the three oracles returned different
verdicts, and returned decompositions whose summands do not actually add
up to the target in the group. A run with a non-zero count is a bug
report, not a measurement, and no metric from it counts.

This is what makes the target safe to optimise against: a faster oracle
that quietly loses decompositions fails the gate.

---

## Baseline, measured

Seed `0x5EED` for structure, `0xB0B` for oracles; 8 targets per instance.

### System structure and first fall degree

16 independent target draws per row. The fall degree is a property of
the system, and the system depends on the target abscissa `x(R)`, so a
single draw is a noisy statistic — `n = 9, m = 2` falls at `D = 2` on 7
draws and shows no fall to `D = 3` on the other 9. Anything built on
this number has to average; `ffd_summary` does.

| n | ℓ | m | vars | eqs | deg | FFD min | FFD max | no fall | mean syz D=2 |
|--:|--:|--:|-----:|----:|----:|--------:|--------:|--------:|-------------:|
| 7 | 3 | 2 | 6 | 7 | 2 | 2 | 3 | 0/16 | 0.38 |
| 7 | 3 | 3 | 16 | 14 | 3 | 2 | 3 | 0/16 | 0.06 |
| 7 | 3 | 4 | 26 | 21 | 3 | 2 | 3 | 0/16 | 0.06 |
| 9 | 6 | 2 | 12 | 9 | 2 | 2 | 2 | 9/16 | 0.44 |
| 9 | 6 | 3 | 27 | 18 | 3 | 3 | 3 | 0/16 | 0.00 |
| 9 | 6 | 4 | 42 | 27 | 3 | 3 | 3 | 0/16 | 0.00 |
| 15 | 4 | 2 | 8 | 15 | 2 | 2 | 3 | 0/16 | 0.94 |
| 15 | 4 | 3 | 27 | 30 | 3 | 3 | 3 | 0/16 | 0.00 |
| 15 | 4 | 4 | 46 | 45 | 3 | 3 | 3 | 0/16 | 0.00 |
| 21 | 6 | 2 | 12 | 21 | 2 | 2 | 3 | 0/16 | 0.50 |
| 21 | 6 | 3 | 39 | 42 | 3 | 3 | 3 | 0/16 | 0.00 |
| 31 | 5 | 2 | 10 | 31 | 2 | 2 | 2 | 0/16 | 10.00 |
| 31 | 5 | 3 | 46 | 62 | 3 | 3 | 3 | 0/16 | 0.00 |
| 63 | 6 | 2 | 12 | 63 | 2 | 2 | 2 | 0/16 | 35.00 |

Three things the table says:

- **No system falls later than `D = 3`,** anywhere on the ladder, at any
  `m` that fits. That is H1.
- **`m = 2` is not the interesting case.** Those systems have `n`
  equations in `2ℓ` unknowns, so they are massively overdetermined by
  `n = 31` (mean 10 syzygies at `D = 2`) and `n = 63` (35). They fall
  immediately because they are nearly linear-algebra problems, not
  because the Semaev structure is weak.
- **Chained systems (`m ≥ 3`) fall at exactly `D = 3` with zero
  syzygies at `D = 2`,** consistently, from `n = 9` up. That is the
  regime the attack actually needs, and it is well-behaved so far.

Comparable numbers for the **full-field** `2n`-variable descent are in
`RESEARCH_FFD_MEASUREMENT.md` (FFD 3 throughout, single draw). The
subspace restriction is what lets `m = 2` fall at 2 — the same
restriction that makes the system solvable at all.

### Oracle cost (8 targets each, `n ≤ 24` so a curve exists)

Medians are split by verdict class; `—` means the class was empty.

| curve | n | ℓ | m | \|F\| | vars | oracle | found | refuted | med found | med refuted |
|:------|--:|--:|--:|------:|-----:|:-------|------:|--------:|----------:|------------:|
| K_0 | 7 | 3 | 2 | 1 | 6 | search | 0 | 8 | — | 0.005 ms |
| K_0 | 7 | 3 | 2 | 1 | 6 | matrix-f4 | 0 | 8 | — | 0.072 ms |
| K_0 | 7 | 3 | 2 | 1 | 6 | sat | 0 | 8 | — | 0.160 ms |
| K_1 | 9 | 6 | 2 | 73 | 12 | search | 8 | 0 | 0.018 ms | — |
| K_1 | 9 | 6 | 2 | 73 | 12 | matrix-f4 | 8 | 0 | 1.642 ms | — |
| K_1 | 9 | 6 | 2 | 73 | 12 | sat | 8 | 0 | 0.731 ms | — |
| K_0 | 9 | 6 | 3 | 55 | 27 | search | 8 | 0 | 0.015 ms | — |
| K_0 | 9 | 6 | 3 | 55 | 27 | matrix-f4 | 8 | 0 | 27.00 ms | — |
| K_0 | 9 | 6 | 3 | 55 | 27 | sat | 8 | 0 | 131.4 ms | — |
| K_1 | 15 | 4 | 3 | 31 | 27 | search | 0 | 8 | — | 3.97 ms |
| K_1 | 15 | 4 | 3 | 31 | 27 | matrix-f4 | 0 | 8 | — | 74.91 ms |
| K_1 | 15 | 4 | 3 | 31 | 27 | sat | 0 | 8 | — | **25 345 ms** |

Search wins everywhere reachable, by two to four orders of magnitude.
That is expected and is not the thing to fix: search costs `|F|^{m−1}`
and cannot reach useful `m` either.  The point of the algebra is that
its cost tracks the system, not `|F|` — which only starts to matter at
an `m` nobody can currently run.

The last three rows are the ones that changed the plan: same variable
count as the rows above them, no decomposition to find, and SAT goes
from 131 ms to 25 seconds.

---

## Pre-registered hypotheses

Each has a falsifier and a bar.  State which you tested and what
happened, in the autolab entry, even when the answer is "no change".

**H1 — the fall degree stays ≤ 3 across the ladder.**
Over 16 draws per instance, every system falls at `D = 2` or `D = 3`,
and every chained (`m ≥ 3`) system falls at exactly 3.
*Falsifier:* any `(n, m)` on the ladder with `fall_min ≥ 4` at
`d_max = 4`, over at least 16 draws.
*Bar:* test the full ladder to `n = 63` at `m ∈ {2, 3, 4}` (those that
fit) and report min/max/no-fall.  Report the distribution, never a
single draw — the spread at `n = 9, m = 2` (7 falls, 9 non-falls) is
what that rule is there for. If H1 holds out to `n = 63`, that is a
publishable-shaped empirical statement about subspace-restricted Semaev
systems; the existing literature measures the full-field case.

**H2 — the 64-variable cap, not the algebra, is what stops the sweep.**
**FALSIFIED, 2026-09-08.** The cap is real (`koblitz_groebner::MAX_VARS`,
a `u64` monomial mask) but it is not what binds. Solve cost saturates at
**27 unknowns**, less than half the cap:

| instance | vars | eqs | eq/var | class | search | F4 | SAT |
|:---------|-----:|----:|-------:|:------|-------:|---:|----:|
| K_0/F_2^9, m=3 | 27 | 18 | 0.67 | found | 0.02 ms | 27 ms | 131 ms |
| K_0/F_2^9, m=3 | 27 | 18 | 0.67 | refuted (raw `x_R`) | — | 6.2 s* | **361 s** |
| K_1/F_2^15, m=3 | 27 | 30 | 1.11 | refuted | 3.97 ms | 75 ms | **25.3 s** |
| K_1/F_2^15, m=2 | 8 | 15 | 1.88 | refuted | 0.21 ms | 0.23 ms | 1.3 ms |

`*` node budget exhausted — F4 did not finish, it gave up.

Widening the mask would let `n = 31, m = 4` (82 unknowns) and
`n = 63, m = 3` (81) be *built*. Nothing suggests they could be
*solved*: they are 3× the variable count at which SAT already needs
minutes. **Do not spend the refactor** — it touches `F2BoolMono`, shared
by ten modules, to buy instances that will not finish.

The primary metric stands, but its justification changes: driving
`unknowns` down matters because solve cost explodes in it, not because
of an arbitrary cap.

**H2′ (replacement) — refutation is the expensive case, and it is where
the oracles differ.**
Finding one root among many is cheap; proving no root exists is not.
At the same 27 unknowns, `K_0/F_2^9` (all targets decompose) costs F4
27 ms, while `K_1/F_2^15` (no target decomposes) costs SAT 25 s — a
340× spread between the two regimes on the same variable count, and a
**340× spread between F4 and SAT on the refutation**.
*Falsifier:* an instance where refutation is not the dominant cost, or
where the two oracles' refutation costs are within 2× of each other.
*Bar:* report `median_found_ms` and `median_refuted_ms` separately —
`OracleRun` now does — for every instance. A single median over both
regimes is a number that describes neither, and it is what hid this
result until 2026-09-08.

**H3 — `S₄` links cut the chaining term roughly in half.**
`S₃` chaining spends one intermediate per extra summand: `m − 2` of
them. A link built on `S₄` absorbs two new summands, needing
`≈ ⌈(m − 3)/2⌉`. At `n = 31, ℓ = 5, m = 7`: 190 unknowns → ≈ 97.
`binary_semaev::binary_semaev_s4` already exists (as a resultant of two
`S₃`s) and needs a symbolic counterpart in `koblitz_groebner`.
*Falsifier:* the `S₄` system disagrees with exhaustive search on any
target, or its higher degree costs more solve time than the variables
saved buy back.
*Bar:* H3 is confirmed if some `(n, m)` that misses the 64-variable
budget under `S₃` fits under `S₄`, solves, and agrees with search.
Note the trade: `S₄` raises the system degree, which raises the
Macaulay column count — measure both, do not assume.

**H4 — F4 beats SAT on chained instances.**
**Supported and sharpened, 2026-09-08.** The gap is concentrated in
refutations: on finds at `n = 9, m = 3` it is 27 ms vs 131 ms (5×), on
refutations at `n = 15, m = 3` it is 75 ms vs 25.3 s (**340×**). The
CDCL solver has no XOR-constraint Gaussian elimination (`semaev_sat`'s
own header says so) and these systems are XOR-dominated, which is the
obvious suspect.
*Falsifier:* SAT wins at any larger `m` or `n`.
*Bar:* three instances at `m ≥ 3` with a consistent winner. If SAT wins
as `m` grows, that inverts which engine deserves the optimisation
effort — worth knowing early, and cheap to test.

**H5 — no algebraic oracle beats exhaustive search anywhere reachable.**
The pessimistic null. Search is 0.004–0.02 ms per target across the
whole benchmarked range; the algebra is 0.07–136 ms.
*Falsifier:* any instance where an algebraic oracle's median is lower,
with `disagreements == 0`.
*Bar:* H5 stands until falsified. Falsifying it is the single most
interesting outcome available from this target, because it would be the
first measured point where the sub-exponential machinery actually pays
on a subfield curve. It most plausibly falls at large `m` with a small
`ℓ` (`n = 63, ℓ = 6` is the best candidate rung) — which is exactly
what H2 and H3 are trying to make reachable.

---

## Ranked next steps

Reordered 2026-09-08 after H2 was falsified. Widening the monomial type
was #1; it is now struck out, because the instances it unlocks cannot be
solved anyway.

1. **Cut refutation cost** (H2′, H4). This is where the time goes and
   where the two oracles differ by 340×. Two concrete moves:
   - **XOR-native propagation in the CDCL solver.** The encoding is
     XOR-dominated and the solver reasons about parity constraints only
     through their CNF expansion. This is the single change most likely
     to move SAT's refutation numbers, and it helps every other
     `semaev_sat` user too.
   - **Make F4's refutations cheaper**, since it already wins: the
     `n = 9, m = 3` raw-target case exhausted the node budget rather
     than returning, so splitting is doing work the algebra should.
2. **Raise the decomposition probability so refutations are rare.**
   A relation search that refutes most targets is paying the expensive
   case almost every time. `|F|^m / m!` against `r` is the knob;
   `subspace_ladder` plus `bench_instance` can map where it sits. This
   may matter more than making refutation faster.
3. **Symbolic `S₄` links** (H3). Still worth doing — it cuts the
   variable count, and after H2 we know variables are expensive for
   real reasons rather than for a cap.
4. **Sparse Macaulay reduction.** `matrix_f4_f2` is dense and capped at
   `MAX_F4_ROWS`/`MAX_F4_COLS`; the `m = 3` systems already reach 6022
   columns at `D = 3`.
5. **Extend the ladder past `n = 63`** for *structure only* — the FFD
   measurement is cheap and `n = 127` (ℓ = 7) is the first rung that
   resembles a deployed curve. Do not expect to solve there.
6. ~~**Widen the monomial type.**~~ Falsified as a priority by H2: it
   buys instances that will not finish. Revisit only if refutation cost
   comes down by orders of magnitude first.

## How to run

```bash
cargo run --release --example koblitz_scaling_bench            # tables
cargo run --release --example koblitz_scaling_bench -- --json  # diff against a baseline
cargo test --release --lib koblitz                             # gates
```

The harness is deterministic: same seeds, same numbers. Report the
primary metric, the three secondaries, which hypotheses were touched,
and the disagreement count — which must be zero.

## What would count as finishing this thread

Any one of:

- H5 falsified: a measured instance where algebra beats search.
- H1 falsified: a fall degree that grows, which would say the
  subspace-restricted systems are not as benign as the first data
  suggests, and would matter to the FFD controversy directly.
- The primary metric brought under 64 at `n = 31` *and* `n = 63`, with
  a solve and a clean gate — that is `m ≈ n/ℓ` reached for the first
  time on a curve of non-trivial size.

None of these threatens a deployed curve, and none is claimed to. The
thread is about establishing where the crossover is, with measurements
rather than extrapolation.

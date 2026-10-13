# Protocol: index calculus on the `F_q` plane quartic against `rho.negation` on `E0(F_{q³})`

Frozen 2026-10-12, before any solver exists or any run. Conductor T-2280.
**This is a plan. Nothing here has been implemented or measured.**
Toy-field measurement on public synthetic instances with planted scalars:
no sub-rho claim is made in advance, none about deployed curves, and the
instances are far below any size of interest. Follows
`research/ghs_cubic_genus3_build_20261012/` (PR #1761, head
`adf5d9b3c99443af74954b6fa0e77c735cba2e68`).

## Question

For `E0` over `K = F_{q³}` with `#E0(K) = 4ℓ` and its smooth plane quartic
`C` over `F_q` (built as in PR #1761):

1. Is the transfer `φ : E0(K) → Jac(C)(F_q)` a homomorphism injective on
   `E0(K)[ℓ]`? This is tested by the separately preregistered checks C1-C4
   below, not by whether a scalar comes back.
2. What does the whole method cost, cold, from the public inputs to the
   returned scalar, against the strongest applicable single-target rho on
   the same point, in one counted unit?

## Workload identities (frozen)

**Panel.** Exactly the eleven records of PR #1761, no others. The panel
contains no field above `q = 11`; the plan does not claim instances that
do not exist. A larger-`q` panel would need its own build and validation
and is a different plan.

| `q` | `u` field indices `(b, c)` | `#E0(K)` | `ℓ` | roots of `u` in `K` |
| --- | --- | --- | --- | --- |
| 5 | (5, 30) | 148 | 37 | 0 |
| 5 | (5, 51) | 148 | 37 | 2 |
| 7 | (7, 52), (7, 54) | 356 | 89 | 0 |
| 7 | (7, 110), (7, 118), (7, 129) | 356 | 89 | 2 |
| 11 | (11, 130) | 1348 | 337 | 0 |
| 11 | (11, 172) | 1324 | 331 | 2 |
| 7 | (49, 8), (49, 9) | 332 | 83 | 0 |

**Field.** `K = F_{q³}` exactly as `Fld::build(q, 3)` in
`research/ghs_cubic_genus3_build_20261012/genus3_build.rs` (SHA-256
`63f7367384ea98643d1acf21f4622450ed5b60ac108bbdcedfc6e6b010371d79`): the
lexicographically first monic cubic modulus over `F_q` (coefficients as
base-`q` digits of the field index) admitting a primitive element, and the
first such primitive element `α` by the tool's search order; element
index `0` is zero and index `e ≥ 1` is `α^(e-1)`; the field index of an
element is the base-`q` integer of its coordinates in the basis
`1, α, α²`. `F_q` is the prime field `Z/q`. `u = t² + b·t + c` with `b, c`
the elements of the stated field indices. The follow-on PR commits an
instance manifest (one JSON record per instance) with the modulus, `α`,
`u`, `u'`, `u''`, the `E0` Weierstrass model `(a4, a6)` and the
map from the quartic model `y² = u·u''` to it (origin: the point at
infinity with `y/t² = +1`), the quartic's coefficients as in #1761's
`quartic_Fq`, and the SHA-256 of the manifest record.

**Subgroup and generator.** `G = [4]P*`, `P*` the first point of
`E0(K)` in the order (field index of `x`, then field index of `y`) with
`[4]P* ≠ O`. Then `G` has order `ℓ` (prime), is deterministic and is
recorded in the manifest.

**Targets.** Per instance, ten workloads `j = 1..10`; `Q_j = k_j·G` with
`k_j` the `j`-th value in `[1, ℓ-1]` of the ecbench workload seed stream
(`seed = 20261012`, stream derivation as the spec records it). The solver
receives `E0`, `C`, `G`, `Q_j` and never `k_j`; ecbench's runner verifies
the returned scalar in its own process.

**Eligible auxiliary points and the F_{q⁶} limit.** PR #1761 shows the
two preimages in `C` of `P = (t0, y0) ∈ E0(K)` are `K`-rational only when
`(u·u')(t0)` is a non-zero square; otherwise they are conjugate over
`F_{q⁶}`, and neither this plan nor ecbench has `F_{q⁶}` arithmetic or a
charge for it. The plan therefore does **not** build a general divisor
path. The index calculus only needs `φ` at auxiliary points `R0 = rG`,
`R1 = R0 + G`, `R2 = R0 + Q_j`; the first `r ≥ 1` (at most 64) with all
three **eligible** is used, and the search for it (K-operations and square
tests) is charged to the workload. The target `Q_j` itself is arbitrary
and need not be eligible. A workload with no eligible `r` is a recorded
failure. Eligibility of ℓ-subgroup points is not known from #1761 (it
tested all of `E0(K)`); the manifest records, per instance, the eligible
fraction of `E0(K)[ℓ]`. An instance where it is zero is dropped and
reported. No claim is made about transfer for ineligible points.

## Methods

**Index calculus on `C` (candidate ids `ic.quartic_full`, `ic.quartic_2lp`;
new ids, never reusing an existing one).**

- Factor base: points of `C(F_q)`. Unknowns: `x_v = λ(P_v − P_0)` and
  `h = λ(H − 4P_0)` for the projection `λ` to the `ℓ`-part, `x_0 = 0`.
- Relations: a line through two base points meets `C` in two more
  points; when both are `F_q`-rational, `x_a + x_b + x_c + x_d = h`.
  The residual pair comes from the restriction of the quartic to the line
  (three evaluations and one square root), then a lookup.
- `ic.quartic_full`: `B = C(F_q)` (all `q + 1` points are materialized);
  dense elimination modulo `ℓ` on `q + 1` unknowns.
- `ic.quartic_2lp`: `B` the first `⌈2.83·√(q+1)⌉` points of `C(F_q)`;
  relations with one or two points outside `B` go through a signed
  union-find (double large prime); cycles give relations among `B`;
  dense elimination on `|B|` unknowns. All pairs of `B` are processed.
- Target: with eligible `R0, R1, R2`, `φ(R)` is an `F_q`-rational
  effective divisor `D_R` of degree 6 over six `K`-points. A cubic through
  `D_R` and three base points meets `C` in three more points; when all are
  `F_q`-rational and resolved, `λ(D_R − 6P_0) = 3h − Σx`. At most 400
  cubics per `R`.
- Logarithm: `y_i = λ(D_{R_i} − 6P_0)`, `k = (y2 − y0)/(y1 − y0) mod ℓ`.

**Reference arm.** `rho.negation`, the repository's strongest applicable
single-target reference on an odd-characteristic curve
(`docs/ecbench/README.md`, section 3 "Boundaries"), on `E0(K)[ℓ]`, same
workload, same seed. No arm "without negation" is preregistered. Its
canonicalization per step is charged as in the pinned conversions below.
The ecbench generic floor `floor_s = √(π/2A)` with `A = 2` is recorded
with every record, as in every other ecbench table.

## Execution vehicle: ecbench, not a standalone executable

New cross-method comparisons run as registered ecbench arms
(AGENTS.md section 12; `docs/ecbench/README.md`): one-unit, one-target,
sealed records, replay-audited. This plan therefore specifies an ecbench
extension (skills `ecbench-extend`, `ecbench-measure`, `ecbench-bounds`),
not a separate `quartic_ic.rs`:

- A new curve construction: a `workload::CurveSpec` variant giving
  `E0(F_{q³})` explicitly (modulus, `(a4, a6)`, `ℓ`, `G`) from the
  manifest, with its group arithmetic over `F_{q³}`; the instance's ICV1
  slug and EC1 alias registered in the same follow-on PR
  (`scripts/build_curve_registry.py`; AGENTS.md section 11). No slug is
  written here because none exists yet.
- The two IC methods as registered methods on `CountedGroup`/`ic_framework`
  plug-ins: deterministic in (workload, seed); they take `Q_j` and the seed
  and never the scalar; set-up and search on separate ledgers (phases);
  every lookup counted in a `*_uncharged` counter unless priced as below;
  a new id for each variant; never altered after a session uses it.
- Spec `ecbench.spec/v1`, one session, `ecbench verify --replay N`, the
  receipt's SHA-256 as the replay certificate; sessions committed under
  `research/ghs_cubic_quartic_ic_20261012/sessions/`.
- Wall time is not a result (see Not measured).

## Boundary and unit, stated before measuring

- **Unit.** `S = total counted cost / √ℓ`, whole method, cold, from the
  public inputs to the returned scalar, with the field multiplications
  behind group operations counted where the curve counts them.
- **Pinned charges (conventions, not measurements).** `F_q`: mult 1, add 1,
  inversion 10, square root `⌈1.5·log₂ q⌉`. `K`: mult 9, add 3, inversion
  90, square root `9·⌈1.5·log₂ q³⌉`, Frobenius `9·⌈1.5·log₂ q⌉`. Modulo
  `ℓ`: mult 9, add 3, inversion 90.
- **Lookups and canonicalization (pinned, both arms).** Each hash or table
  lookup or insertion, in either arm (relation lookup, union-find step, rho
  visited-point table), is charged 4 `F_q`-multiplication equivalents; each
  `rho.negation` canonicalization is charged as one `K` addition and one
  comparison (3 + 1). These conversion costs are preregistered here, before
  any measurement, so that neither arm is a bare lower bound.
- **If the harness cannot carry the pinned charges** (ecbench sets
  `cost.lower_bound` and `bounded: true` whenever anything is unpriced,
  `docs/ecbench/README.md` section 3), then rows from both arms are
  counted-work diagnostics. A ratio of two lower bounds is not a bound on
  the ratio of full costs; such a row cannot earn "below rho".
- **Setup-phase prediction (not a floor).** For `ic.quartic_full`,
  setup costs about `(q+1)·c_pt` with `c_pt` a future measured per-point
  cost. This is a prediction of one phase: its coefficient is measured
  later and it binds only if every point is materialized. It is labelled
  a prediction everywhere.
- **Floors (derived, per variant).** A linear system needs at least as
  many independent equations as unknowns, and each relation costs at
  least one `F_q` operation (one evaluation of the restricted quartic):
  `S_floor(ic.quartic_full) = (q+1)/√ℓ` and
  `S_floor(ic.quartic_2lp) = ⌈2.83·√(q+1)⌉/√ℓ`. These are deliberately
  loose counting bounds, far below any plausible cost, and say nothing
  about whether the method crosses rho. The ratio `S / S_floor` is
  reported beside `S`.

## Pre-registered predictions

- **P1 (pipeline correctness).** Every index-calculus run that returns,
  returns `k`. A wrong scalar invalidates the composite pipeline only:
  a defect in divisor arithmetic, relation generation, union-find,
  elimination or target descent gives the same observation, so a wrong
  scalar alone does not falsify the transfer derivation of PR #1761.
  The transfer is judged by the checks below.
- **P2 (`ic.quartic_full`).** `S_ic.quartic_full / S_rho` increases with
  `q` across `q ∈ {5, 7, 11}` (dense elimination on `q+1` unknowns is
  cubic in `q`; rho is `q^{3/2}`). Three values of `q` are a weak test,
  and the panel supports no trend claim beyond saying so.
- **P3 (`ic.quartic_2lp`).** No prediction on the ratio's trend or on a
  crossing.

**Localizing checks (separately preregistered, run before any solve is
timed or compared).**
C1 `#Jac(C)(F_q)·D ~ 0` for random degree-0 `D`, and associativity of the
class arithmetic, on all eleven quartics.
C2 `φ(P+Q) ~ φ(P) + φ(Q)` for random eligible `P, Q`.
C3 `ℓ·φ(P) ~ 0` and `φ(P) ≁ 0` for `P` a generator of `E0(K)[ℓ]`, at least
one case per instance where eligible points of the `ℓ`-subgroup exist.
C4 each emitted relation re-verified by the independent linear-equivalence
routine on a fixed sample. A failing check localizes the failure to that
component; P1 alone never does.

## Decision rule (per paired one-target workload)

The unit of verdict is one paired workload: one instance, one target
`Q_j`, one IC arm against `rho.negation` on the same point, same seed and
resource envelope (AGENTS.md, "Primary comparison uses one target"). Each
workload has its own result row: returned scalar verified or not, `S`
for both arms, the ratio `S_ic / S_rho`, `S / S_floor`, the `bounded`
flag, and the failure reason if any. Failures (rank not reached, no
decomposition within 400 cubics, no eligible `r`, wrong scalar) are
retained in the panel, not dropped.

- A workload's IC row is "below rho on this workload" only if its answer
  is verified, both arms are fully charged (`bounded: false`), and
  `S_ic < S_rho` with the paired bootstrap interval the harness records
  excluding 1.
- The report states the number of workloads with that outcome out of the
  total, per variant and per instance, and the table of every row. It
  makes no pooled "below rho on these instances" claim, and no crossover
  claim, from three values of `q`. Any cross-target mean or other batch
  summary is secondary and descriptive only; it never replaces or decides
  the per-target rows.
- Rows with `bounded: true` stay counted-work diagnostics and cannot earn
  "below rho".

## Not measured

Wall time; memory; sparse linear algebra; `F_{q⁶}` arithmetic and
ineligible auxiliary points; any field above `q = 11`; the hyperelliptic
`m = 3` covers; `n ≠ 3`; a multi-target DLP run.

## Language and report path

Rust only (AGENTS.md, "Implementation language: no Python"), inside the
repository's library and ecbench; no Python for any execution, harness,
checker, or result generation. Graphs are SVG written by a native Rust
renderer from the recorded session values (ecbench records only, nothing
recomputed), and the report PDF is built from a Typst source with the SVG
embedded, per the `research-visuals` skill; the PNG render, where the
benchmarks rule requires one, comes from the same SVG with a native
rasterizer. The report records the command that produced each file.

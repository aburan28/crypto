# Lead 5b: the `D₄`-symmetrised `S₅` as four weighted unknowns

2026-10-15, lane `d4-s5`.

- Protocol: [`PROTOCOL.md`](PROTOCOL.md), committed in `87fe78b63` before any
  measured run.
- Code: `src/cryptanalysis/d4_weighted_s5.rs`, `examples/d4_weighted_s5.rs`,
  `examples/d4_weighted_s5_summary.rs`.
- Report: [`REPORT.pdf`](REPORT.pdf).
- Run window: 2026-10-15 06:38Z–07:43Z (`raw/started_utc.txt`,
  `raw/finished_utc.txt`).
- Host: Claude Code cloud container, 4 vCPU Xeon, shared by five lanes.
  Each run used one thread, `nice 10` and a 2.5 GB address-space cap.
- Class: **engineering** (a constant of one phase). No end-to-end `S` was
  measured, and none is claimed.

## 0. The lead was not "never run"

The encyclopedia's §8 item 5 says the `D₄`-symmetrised `S₅` was "never run".
That entry is superseded by PR #1785, which ran it.

- It copied candidate C8 of `research/pdp_root_finder_requirements_20261012`.
- C8 was built and run on 2026-10-13: PR #1785, commit `0846e6a12`,
  `research/gaudry_quartic_d4_20261013/`.
- C8 measured `C₄ = 2.54·10¹⁰` at `p = 269`, which is `48×` below the plain
  form and not the `≈ 630×` projected.

C8 wrote the system as five unweighted unknowns `(σ₁, σ₂, σ₃, σ₄, π)` with
`π² = σ₄`. This study builds the presentation the GLV note §9.4 named and
nobody built: a **weighted Macaulay solve** in four unknowns. It also adds
what C8 left open:
- several field sizes;
- brute-force ground truth;
- the lift-back cost, charged;
- the handover at the cofactor E18 measured.

## 1. Construction (derived; see PROTOCOL §1)

**Curve.** The curve is over `F_{p⁴}` with a rational `T = (x₀, 0)` and
`f'(x₀) = c²`, `c ∈ F_p^*`.

**Coordinate and base.** FGHR's coordinate is
`Y = (x − x₀ − c)/(x − x₀ + c)`. It satisfies `Y(P + T) = −Y(P)` and
`Y(−P) = Y(P)`. The base is `{P : Y(P) ∈ F_p \ {±1}}`, which is E18's
`τ_T` line.

**Group.** `D₄ = (Z/2)³ ⋊ S₄`, of order 192, acts on the four summands'
`Y_i`. It permutes them and changes an even number of signs, which is the
same as translating an even number of summands by `T`.

**Invariants.** They are `σ_k = e_k(Y²)` for `k = 1, 2, 3` and
`π = Y₁Y₂Y₃Y₄`. `H = S₅ · Π(1 − Y)⁸` lies on the 105 monomials
`σ₁^a σ₂^b σ₃^c π^d` with `2(a+b+c) + d ≤ 8`.

**System.** Four Weil components over `F_p` in `(σ₁, σ₂, σ₃, π)`, with
weights `(2, 2, 2, 1)` and **weighted degree 8**.

**Bézout `4096 → 512`, derived two ways (labelled derived):**
- weighted Bézout: `8⁴/(2·2·2·1) = 512`;
- orbit count: the `S₄` form's 4,096 solutions are `S₄`-orbits, and
  `|D₄|/|S₄| = 8` of them make one `D₄`-orbit, so `4096/8 = 512`.

The weighted Macaulay bound is `Σdᵢ − Σwⱼ + 1 = 26`. At that degree the
matrix has **4,200 columns**, against C8's 11,628 and the plain form's
40,920.

**Measured confirmation.** The quotient dimension is 512 on every solved
target at `p ≥ 29`, and on 3 of 5 at `p = 17` (§3).

## 2. Controls (measured, `raw/controls.jsonl`)

- **Negative, the invariance detector.** With the true `c`, `H` lies in
  the `D₄`-invariant span. With a wrong `c`, the interpolation fails at the
  first check point. Both hold at `p = 13` and `p = 269`.
- **Positive, the solver.** Nine planted generic weighted systems
  (`p = 269, 1033, 4093`, 3 seeds each) all gave a `3,979 × 4,200` matrix,
  quotient 512 and the planted root. Each cost
  **`1.386–1.451·10⁹`** multiplications.
- **Positive, the pipeline.** Every planted target that closed returned its
  planted quadruple: 41 of 42 for arm w, 23 of 24 for C8. The exception is
  the degenerate target below.

## 3. Results (measured unless marked; `summary.txt`, `summary.json`)

Units: `F_p` multiplications per solve, counted by `gaudry_quartic`'s
convention. That is the unit of `C₄ = 1.213·10¹²`. Arm w's figures include
any reruns at degree 27.

| p | seed | targets | `C_w` median (range) | `C_c8` median (range) | paired `C_c8/C_w` median (range) | `1.213·10¹²/C_w` median (range) | lift-back / `C_w`, max | w = ground truth |
|--:|--:|--:|:--|:--|:--|:--|--:|:--|
| 13 | 1 | 6 | `1.21·10⁹` (`1.17–2.82·10⁹`) | `2.18·10¹⁰` | 18.0 (7.7–18.6) | 1002 (430–1042) | `4·10⁻⁵` | brute 6/6, MITM 6/6 |
| 17 | 1 | 6 | `1.21·10⁹` (`1.06–2.54·10⁹`) | `2.05·10¹⁰` | 18.5 (8.1–19.3) | 1005 (477–1146) | `7·10⁻⁵` | brute **5/6**, MITM 5/6 |
| 29 | 1 | 6 | `1.33·10⁹` (`1.27–3.10·10⁹`) | `2.39·10¹⁰` | 18.0 (7.7–18.7) | 912 (392–953) | `9·10⁻⁵` | brute 6/6, MITM 6/6 |
| 37 | 1 | 6 | `1.34·10⁹` (`1.29–1.38·10⁹`) | `2.42·10¹⁰` | 17.9 (17.6–18.7) | 905 (881–942) | `9·10⁻⁵` | brute 6/6, MITM 6/6 |
| **269** | 1 | 8 | **`1.411·10⁹`** (`1.36–1.57·10⁹`) | `2.542·10¹⁰` | **18.0** (16.3–18.7) | **860** (772–893) | `1.3·10⁻⁴` | MITM 8/8 |
| 269 | 2 | 12 | `1.385·10⁹` (`1.36–3.33·10⁹`) | — | — | 876 (365–894) | `1.4·10⁻⁴` | MITM 12/12 |
| **1033** | 1 | 8 | **`1.406·10⁹`** (`1.37–1.47·10⁹`) | `2.554·10¹⁰` | **18.2** (17.4–18.7) | **863** (823–888) | `1.6·10⁻⁴` | MITM 8/8 |
| 1033 | 2 | 12 | `1.420·10⁹` (`1.39–1.53·10⁹`) | — | — | 854 (794–871) | `1.7·10⁻⁴` | MITM 12/12 |
| **4093** | 1 | 8 | **`1.437·10⁹`** (`1.37–1.56·10⁹`) | `2.561·10¹⁰` | **17.8** (16.5–18.6) | **844** (778–885) | `2.0·10⁻⁴` | MITM 8/8 |
| 4093 | 2 | 12 | `1.424·10⁹` (`1.37–1.51·10⁹`) | — | — | 852 (806–885) | `2.0·10⁻⁴` | MITM 12/12 |

**Matrix and phases.**
- Arm w's matrix is `3,979 × 4,200` (33 MB); C8's is `13,412 × 11,628`.
- Elimination falls from `1.90·10¹⁰` to `4.8·10⁸`, `÷ 40`.
- Normal forms fall from `6.2·10⁹` to `7.5·10⁸`, `÷ 8`.
- The characteristic polynomial is the same `1.4·10⁸` in both arms
  (`≈ 1.05·512³`), because both have `D = 512`.
- In arm w the normal forms, not the elimination, are now the largest phase.

**Degree-27 reruns.** 5 of 101 solved targets had a degenerate `M_π`
eigenspace and reran at degree 27, at 4,760 columns. All 5 are charged, at
`2.5–3.3·10⁹` each. Four were at `p ≤ 29`, where `π`-values collide often;
one was at `p = 269`, seed 2.

**Branches** (per target, `raw/*.jsonl`):
- The rational invariant solutions are mostly spurious. At the measured
  sizes, 126 of them gave 34 that split into `Y_i²` and 30 rooted.
- The 30 rooted solutions gave 240 candidate sign vectors (8 each) and
  **240 verified decompositions**.
- Every planted target at `p ≥ 29` decomposes into exactly 8 sorted
  quadruples, which are the 8 even `T`-translates. At `p = 13` it is 2, and
  at `p = 17` it is 2–8, because translates leave the base or collide.
- C8 returned the same set as arm w on all 47 paired solved targets.

**Lift-back** (charged: split, `F_p` square roots, `x(Y)`, lifts in
`F_{p⁴}`, and group checks of all 8 translates):
- at most `2.8·10⁵` multiplications;
- at most `2.0·10⁻⁴` of `C_w`.

**The one failure.** `p = 17`, planted target 4. Its quotient has dimension
514, not 512: a degenerate system of the same kind as C8's dimension-517
target. Arm w returned *unsolved* (border unreachable), not a wrong set.
C8 also failed on it. Brute force finds 4 decompositions.

## 4. Against the registered predictions and decision rule

- **W1 (dimension 512, 4,200 columns): held at `p ≥ 29`, not at `p = 17`.**
  - At `p ≥ 29`, every target had dimension 512 and 4,200 columns.
  - At `p = 17`, dimensions were 510, 511, 511, 512, 512 and the unsolved
    514.
  - At `p = 13`, all six had dimension 512.
- **W2 (exactness on every target): not met as registered.**
  - Of 102 targets, 101 closed, and on every one of them arm w's
    decompositions equal brute force (tiny fields) and meet in the middle.
  - The `p = 17` target above is a miss: unsolved and flagged, shared with
    C8, not a wrong answer.
  - PROTOCOL §6 says any disagreement voids arm w's cost numbers. Read
    literally, that rule applies here. So the decision below is reported
    **conditional on that caveat**, and the cost numbers are stated as
    measured, not admitted as a clean verdict.
  - In a pipeline the failure costs one discarded solve. That is one in 102
    here, and C8 saw one in 20.
- **W3 (flat in `p` within `±25 %`): held.** The seed-1 medians are
  `1.411`, `1.406` and `1.437·10⁹`, a spread of 2.2 %.
- **W4 (projection 630×, `C ≈ 1.9·10⁹`):**
  - The measured `C_w` is `1.41–1.44·10⁹`.
  - The generic control is `1.39–1.45·10⁹`.
  - Both sit 25 % *below* the `D^{3.1}` projection.
- **Decision on lead 5b: "projection met".**
  - The median `R_plain` is 860, 863 and 844 at `p = 269`, 1033 and 4093,
    all inside `[315, 1260]`.
  - **Not "accounting":** the lift-back is `≤ 2·10⁻⁴` of `C_w`, under the
    10 % threshold.
  - This verdict carries the W2 caveat above.
  - The plain arm itself was not re-run. It is quoted from §11.17, because
    its 3.42 GB matrix exceeds the cap.
- **Against C8: "advances on C8".** The paired `C_c8/C_w` is at least 2 on
  every one of the 47 paired targets: minimum 7.7 (degree-27 reruns),
  median 17.8–18.5.

So the gap between C8's `48×` and the projected `630×` was the presentation,
not the symmetry. Adjoining `σ₄ = π²` as a fifth unknown keeps the
unweighted grading, and that triples the columns. With weights, the
`D₄` system costs what a generic system of quotient dimension 512 costs.

## 5. The `k = 4` handover (derived, extrapolation; not a measurement)

The frame is E18's `τ_T` row (`RESEARCH_GLV_INVARIANT_FACTOR_BASES.md`
§9.3). The `D₄` system needs rational `T`, so the curves are E18's 2-torsion
curves. Its inputs are measured:
- `h = 4`, which is what E18 found on all 8 curves;
- `y = 106` residuals per relation;
- `r = 0.341` against the negation walk.

`n*` grows as `C⁴`, so `n* = 157.9 − 4 log₂(1.213·10¹²/C)` bits.

| solve | `C` (median, measured sizes) | `n*` |
|:--|--:|--:|
| plain `S₄` form (quoted) | `1.213·10¹²` | `2^{157.9}` (E18) |
| C8, five unknowns | `2.555·10¹⁰` | `2^{135.6}` |
| **w, four weighted unknowns** | **`1.420·10⁹`** | **`2^{118.9}`** |
| w + lift-back | `1.420·10⁹` | `2^{118.9}` |

Reference handovers, all extrapolated:
- plain prime-order method: `2^{151.1}` against the unfolded walk (§11.19),
  `2^{156.1}` against the negation walk (E18);
- E18's `ψ`-folded `j = 0` arm: `2^{146.0}`.

**Reading.**
- On these curves the solve moves the `τ_T` handover **39 bits** earlier
  than the plain `S₄` form (`4 log₂ 854`). That is 17 bits earlier than C8.
- It is also 27–37 bits earlier than the `ψ`-folded and plain prime-order
  figures. That comparison crosses curve families: the 2-torsion curves pay
  `h = 4` and a 4.7× poorer yield, which are already inside the `2^{157.9}`
  baseline.
- The note's `≈ 37` bits for a 630× cut is consistent with this.
- Past `n*`, the method tends to `r = 0.341×` the negation walk on its own
  curves. `r∞` does not depend on the solve, so it is unchanged.

Every `n*` above rests on:
- the exponents `n^{1/4}` for relations and `n^{1/2}` for the linear algebra
  and rho;
- a per-solve `C` measured at `p ≤ 4093`. It is flat there, and its
  `p`-dependence (root finding) is about `10⁻⁴` of it.

No `k = 4` run below rho exists.

C8's own model (rate `1/8`, `h = 2`) put the handover at `2^{135}` with
`C = 2.54·10¹⁰`. E18 found `h = 4` on every 2-torsion curve it built, so the
`h = 2` row of C8's README is probably unattainable for this coordinate.
That is a measured tendency, not a proof.

## 6. What this does not show

- The plain `S₄` arm was not run on these curves or targets: its matrix is
  3.42 GB, over the cap. `R_plain` uses §11.17's `C₄`, which was measured on
  generic, non-2-torsion curves on the `x`-line. §11.17 showed that number
  equals a generic solve of four octics, but it is not a same-target pairing.
- No relation collection, linear algebra or rho was run. The handover is
  derived.
- The `p = 17` degenerate target's cause is not established, nor is that of
  C8's dimension-517 target.
- The module's unit tests run once the coordinator adds the missing fixtures
  to the sparse checkout. `cargo test --release --lib -- d4_weighted
  gaudry_quartic` then passes, 11 of 11.
  - The first version of the small-field test planted `P₀ + P₁ + P₂ − P₃`
    at `p = 13` (instance seed 2). On that instance this sum is `O`, so the
    test panicked before any solve ran.
  - The first nonzero planted sum on that instance is a degenerate system,
    quotient dimension 514. Neither this solve nor C8's closes it, and brute
    force finds one decomposition. This is the same failure as `p = 17`,
    target 4, in §3.
  - It is not a solver defect. The test now uses `p = 29`, instance seed 1,
    where every registered target closed. The measured solve path is
    unchanged (commit after `87de2fa83`).
- No `identity.certificate/v1` was produced for the Bézout count. It is
  derived in PROTOCOL §1 and confirmed by the measured quotient dimensions.

## Cross-check after merging `main`

The shared `cargo-target` could in principle serve a stale library, so the
example was rebuilt in a private target directory. That build is after the
merge of `origin/main`, which changed no solver file this study uses.
Targets 0 and 1 of seed 1 were re-run at `p = 269`, `1033` and `4093`
(`raw/xcheck_own_target_dir.jsonl`). For all 6, `C_w`, the lift-back count
and the decomposition count equal the recorded values exactly.

## Related study on `main`

`research/ic_improve_20261015/04-d4-2lp-fix` (benchmark
`20261015-ic-improve2-04-d4-2lp-fix`) improves the relation stop and linear
algebra of the `D₄` `k = 4` pipeline by a constant, `0.979×`. It does not
change the per-decomposition solve this study measures. Its handover
figures (`2^{116.5}`, `2^{128.9}`) are for the 2LP switch in its own model,
not for `C₄`, so they are not comparable with §5.

## Figures

- [`figures/c_per_solve.svg`](figures/c_per_solve.svg): every target's `C`
  for both arms, the quoted plain level and the projection.
- [`figures/phases.svg`](figures/phases.svg): elimination, normal forms,
  characteristic polynomial and eigenvectors for each arm.
- [`figures/handover.svg`](figures/handover.svg): `n*` in E18's frame.
- Source: [`figures/plot.py`](figures/plot.py), which plots only, from
  `raw/*.jsonl` and `summary.json`.

## Reproduce

The measured binary was built at `87fe78b63`. The later commit only replaces `q % 2 == 0` with `q.is_multiple_of(2)` (a clippy fix with the same behaviour) and adds the summariser. Build:

```sh
cargo build --release --example d4_weighted_s5 --example d4_weighted_s5_summary
# the registered cells, in order (raw/ is append-only)
sh research/d4_symmetrised_s5_20261015/raw/run_all.sh   # paths inside point at this lane's checkout
target/release/examples/d4_weighted_s5_summary research/d4_symmetrised_s5_20261015/raw > summary.txt
python figures/plot.py
```

Conductor's control plane was unreachable (connection refused): no task card.

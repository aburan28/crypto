# Protocol: a cheaper `k = 4` decomposition over `F_{p⁴}` (the `D₄` oracle's `C₄`)

Registered 2026-10-14, before the ladder and the ecbench session below.
Code frozen at commit `2b49f38d1` (branch
`claude/k4-d4-decomp-cost-20261014`, on top of PR #1792's head `96cd7314d`,
which is not merged). Binaries built from it:
`target/release/ecbench` sha256 `67b8cc95…3e22`,
`target/release/examples/gaudry_quartic_d4` sha256 `c818a579…f524`.
Class: **engineering** (a constant in `C₄`; the exponent `q^{2−2/n}` of the
method is unchanged).

## Question

`research/gaudry_quartic_d4_20261013` measured `C₄ = 2.54·10¹⁰` `F_p`
multiplications for one decomposition with FGHR's `D₄` system at `p = 269`
(Macaulay degree 14 in five unknowns `(σ₁, σ₂, σ₃, σ₄, π)` with `π² − σ₄`;
elimination 75 %, normal forms 25 %), and
`research/gaudry_quartic_d4_ecbench_20261013` ran the method end to end at
`p = 29` (`S = 6.49·10⁸` with `fghr-y-line:zero=1`, 574 trials). How far can
the counted cost of one decomposition, and with it `S`, be cut while every
decomposition stays exact, and where does §11.16's crossover with rho move?

## What changed (oracle options, all off by default)

`gaudry-d4:…` gains four switches; with all four off it is the committed
oracle, operation for operation (its counts on the committed `p = 269`
planted solves were reproduced exactly before this protocol; the unit tests
pin the rest).

1. `weighted=1`: substitute `σ₄ = π²` and solve the four equations in
   `(σ₁, σ₂, σ₃, π)` in the weights `(2, 2, 2, 1)` (each `σ_k` has degree 2
   in every `Y_i`, `π` degree 1). The equations have weighted degree 8; if
   their top parts are regular the Hilbert series is
   `(1 + t² + t⁴ + t⁶)³(1 + t + … + t⁷)`, so the weighted Macaulay degree is
   26: a **3,979 × 4,200** matrix against **13,412 × 11,628**. The
   multiplication matrix is `π`'s (weight 1, so its border stays in the
   matrix).
2. `chains=1`: only 44 of the 512 products `π·b` leave the standard
   monomials, so `M_π` is block-companion along the chains `h, hπ, …`. The
   step computes only the 44 border normal forms, the characteristic
   polynomial as `det A(λ)` of a 44 × 44 matrix (every point of `F_p` when
   `p ≤ 513`, otherwise 513 points interpolated and its roots found), and
   the eigenvectors as kernels of `A(λ)`. A degenerate eigenspace is split
   by `σ₁, σ₂, σ₃` on the standard monomials whose products stay standard.
3. `sparsest=1`: each elimination column takes the sparsest candidate row
   as pivot instead of the lowest-indexed.
4. `replay=1`: the first solved system's shape (F4Remake's idea) is replayed:
   rows that reduced to zero are not eliminated. If the replayed rows fall
   short of the recorded rank, held rows are reduced in until it is reached
   (top-up). The recorded rank is the generic one and a specialised system
   cannot exceed it, so the result is the full row space. A mismatch of
   leads or row count runs the full solve, and both are charged.

`charge_unsolved=1` is used on both IC arms below, so every attempt is in
`S`.

## Pilot (before this protocol, not cited as results)

Scratch runs of `examples/gaudry_quartic_d4.rs` chose these options and the
predictions below:
- `p = 269`, planted, 2–4 targets per configuration: committed
  `2.55·10¹⁰`; weighted `1.49·10⁹`; weighted + chains `5.09·10⁸`;
  + sparsest `1.41·10⁸`; + replay `6.9·10⁷`; plain + sparsest `1.0·10⁹`.
- Weighted degrees 22–25 do not close (dimension 502–511); 26 does.
- `p = 613, 1009`: interpolated path, `8.2·10⁷` replayed.
- `p = 29`, 300 uniform targets on the ecbench curve: replay without top-up
  fell back on 42 % of systems; with top-up, mean `6.23·10⁷`; 3 of 300
  systems did not close (dimension 509–511), all three without
  decompositions; the committed oracle closed all 300; the decomposition
  sets were identical on the 297 both closed.
- A two-arm scratch plumbing session in ecbench (20 trials).

## Part A — decomposition ladder (`examples/gaudry_quartic_d4.rs`)

Every solve is checked against `mitm_y_decompositions` (exact enumeration
over the same base). "New" is `--weighted --chains --sparsest --replay`
(`--macaulay 26`); "committed" is the default (`--macaulay 14`). Outputs go to
`raw/`.

| p | curve | new | committed (same targets) |
|:--|:--|:--|:--|
| 29 | `--find 16` (the ecbench curve) | 20 planted, 300 uniform | 20 planted, 50 uniform |
| 269 | `--seed 1` | 20 planted, 50 uniform | 10 planted |
| 1009 | `--seed 1` | 20 planted, 50 uniform | 5 planted |
| 4049 | `--seed 1` | 10 planted, 20 uniform | 3 planted |

Ablation at `p = 269`, `--seed 1`, 10 planted targets each: committed;
`--sparsest`; `--sparsest --chains`; `--sparsest --chains --replay`;
`--weighted`; `--weighted --sparsest`; `--weighted --chains`;
`--weighted --chains --sparsest`; `--weighted --chains --sparsest --replay`.

Statistics: mean `C₄` per configuration with a 95 % bootstrap interval over
targets (10,000 resamples, seed 1). With replay, the first (recorded) solve is
reported separately and the mean is over the replayed ones. Ratios
committed/new use a paired bootstrap over shared planted targets.

## Part B — ecbench session (`specs/SPEC-k4-p29-d4cost.json`)

`fp4_torsion(29, 1, 16)`, two targets (`target_seed` of `SPEC-k4-p29-zero`;
the first workload is that session's, id `Wc87e324798ff`). Arms:
`rho.negation` (reference), `ic.pipeline` + `fghr-y-line:zero=1` with
`gaudry-d4:macaulay=14,charge_unsolved=1` (committed oracle) and with
`gaudry-d4:macaulay=26,weighted=1,chains=1,sparsest=1,replay=1,charge_unsolved=1`
(new), `rho.negation` A/A. One round, `max_trials = 20,000`. Then
`ecbench verify --replay-all`, `ecbench compare` (new vs committed, new vs
rho) and `ecbench table`. `S` is a count (`priced_by = counted`, 97 `F_p`
multiplications per addition), so isolation grades only the wall time; the
ladder may run beside the session, and the grades are reported as recorded.

## Part C — crossover

§11.16's model as used in `../gaudry_quartic_d4_20261013`:
`S/rho = √h · (181 · (p/4) · C₄ / (126 p²) + r∞/4)`, `r∞ = 0.34·(1/4)`
at `φ = 0.73` (0.63/4 at `φ = 1`), with `C₄` as measured at large `p` (the
interpolated path; its root-finding grows like `log p` and is charged at the
measured `p` and extrapolated linearly in `log p`). `n* = p⁴` where the
ratio crosses 1, for `h = 1, 2, 4`.

## Predictions

- **A1 (exact).** Every decomposition set from the new oracle equals the
  enumeration, on every target it closes. On every target both oracles
  close, they return the same rational solutions and decompositions.
- **A2 (`C₄`).** New, replayed mean, in `[4·10⁷, 1.2·10⁸]` at each of
  `p = 29, 269, 1009, 4049`. Committed/new at least 200× at each `p` with
  paired targets.
- **A3 (degenerate).** The new oracle leaves at most 3 % of uniform
  targets unsolved at each `p`.
- **A4 (ablation).** `--sparsest` alone cuts the committed `C₄` at least
  10×. At equal linear algebra (`sparsest + chains + replay`), the weighted
  form is at least 5× below the unweighted one.
- **B1.** Every arm verifies `[d]G = Q`; the audit replays every run
  identically.
- **B2.** The two IC arms' trial counts agree within 5 % per target (both
  oracles are exact on the same walk; they can differ only through
  unsolved systems and the witness a solve returns).
- **B3.** `S_committed / S_new` lies in `[250, 500]` on each target; on
  target 0, `S_new` is about `1.7·10⁶` (against `6.49·10⁸` sealed).
- **C1.** `n*` at `h = 2` falls in `[2^{95}, 2^{110}]` (from `2^{135}`).

## Decision and stop rules

A2 and B3 holding is a constant-factor result on `C₄` and on `S`, reported
as such; the crossover is a model, not a measurement. A miss on A1 at any
target stops the claim: the cause is found and reported, and no cost from
that configuration is quoted. A miss on A2–A4, B2–B3 or C1 is reported as
a miss with the measured value. No exponent claim is made.

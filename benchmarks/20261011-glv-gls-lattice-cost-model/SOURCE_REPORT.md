# Isogeny-loop GLV, round 2: open leads, an erratum, and 4-dimensional GLV-GLS

**Date.** 2026-10-11.  **Question.** After round 1 (`RESEARCH_ISOGENY_LOOP_GLV.md`,
PR #1644, evidence in `../isogeny_loop_glv_20261010/`), which leads around
Koshelev–Sanso isogeny-loop endomorphisms (ePrint 2024/1985) are still open for
GLV/GLS-type scalar multiplication on Koblitz, binary and prime curves?
**PDF.** [`REPORT.pdf`](REPORT.pdf) (source [`REPORT.typ`](REPORT.typ)), with the
figures in [`figures/`](figures/) (SVG, and their plotting script).
**Run reports.** [`reports/run-20261011T2258Z-degree-bound/`](reports/run-20261011T2258Z-degree-bound/),
[`reports/run-20261011T2258Z-glvgls-lattice/`](reports/run-20261011T2258Z-glvgls-lattice/).

Status words used below: **measured** (a native run recorded it), **derived**
(follows from stated facts by the argument given), **proposed** (not tested).

## 0. Two corrections to round 1, stated first

1. **Language rule.** Round 1's sweeper (`scripts/isogeny_loop_sweeper.py`) is
   Python with cypari2. AGENTS.md ("Implementation language: no Python") had
   forbidden Python for research result generation since 2026-10-03, so
   PR #1644 breached it. The round-1 outputs stay as frozen evidence with
   Python provenance. All new computation in this round is native Rust:
   `examples/isogeny_loop_degree_bound.rs` and `examples/isogeny_loop_glv_gls_lattice.rs`.
   The sweeper is now marked legacy. Porting it needs a native
   subexponential class-group algorithm for 130-bit discriminants, which the
   crate does not have (`src/isogeny/class_group.rs` enumerates forms up to
   |D| of about 10^6).
2. **Erratum: Verschiebung pricing in characteristic 2.** Round 1 charged both
   directions of a 2-step on a binary curve 2 M. Only the Frobenius direction
   F (two squarings) is that cheap. The conjugate direction is the
   Verschiebung V = [2] ∘ F⁻¹, a doubling composed with square roots.
   - *Original statement* (`RESEARCH_ISOGENY_LOOP_GLV.md` §4.3): the
     folklore split λ = t − 1 via Vⁿ has ratio 0.49 under the 2 M/step model
     and is "not established as a gain".
   - *Correction* (derived): Vⁿ = [2ⁿ] ∘ π⁻¹, and π acts as 1 on E(F_q), so Vⁿ
     is n doublings plus square roots. It is spent to save n/2 doublings, for a
     ratio ≥ 2(c_dbl + 2)/c_dbl = 2.5 at c_dbl = 8. That is a provable loss.
   - *Effect on the other round-1 conclusions*: none. Raising a step weight
     can only raise the cost of every loop, so every "no loop under budget"
     and every "cheapest loop ≥ x" result of round 1 still holds a fortiori.
     The round-1 binary "cheapest loop" costs are lower bounds on the corrected
     costs.
   - The sweeper's cost function is corrected in this PR so the legacy tool no
     longer reproduces the error. Its Python cross-check runs of this
     correction are not used as evidence.

## 1. Lead: a class-group-free lower bound (derived; native table; measured inputs)

A loop of prime-degree steps l_1 … l_m is an endomorphism α of degree
N(α) = ∏ l_i. Its cost is Σ c(l_i) ≥ κ · log₂ N(α), where κ = min_l c(l)/log₂ l.
A non-scalar α = (A + B√D)/2 has B ≠ 0, so N(α) ≥ d_min(D) = ⌊(|D| + 1)/4⌋,
and d_min only grows in sub-orders. Therefore every loop endomorphism of a curve
whose Frobenius has fundamental discriminant D_K costs at least
**κ · log₂ d_min(D_K)**, with no class group needed.

| per-step model | κ (M per bit of degree) | source |
|:--|--:|:--|
| projective Vélu, 7.5 l (the paper's §2.2) | 14.196 (l = 3) | Koshelev–Sanso |
| affine Kohel, 2 l + 2, inversions free | 5.047 (l = 3) | round 1 |
| optimistic x-only, 3-isogeny eval 4M + 2S counted S = M | 3.786 (l = 3) | SIKE-style formulas |

Measured inputs, from `degree_bound.md` / `.json`: for every curve the tool
re-verifies t² − 4q = f² · D_K in big integers against the frozen round-1 records.

| curve family | log₂ d_min ≈ | bound / saving, projective | affine | optimistic x-only |
|:--|:--|--:|--:|--:|
| random prime and binary challenge curves (79–131-bit r) | log₂ r | 3.27–3.51 | 1.16–1.25 | 0.87–0.94 |
| CryptoPro-A (256-bit, class group never computed) | 255.9 | **3.55** | **1.26** | 0.95 |
| CryptoPro-B (D = −619) | 7.3 | 0.10 | 0.04 | 0.03 |

**Reading (derived).** For a random curve |D_K| ≈ q ≈ r, so the ratio is about
2κ/c_dbl. A profitable loop on a random curve therefore needs a per-step
formula cheaper than c_dbl/2 = 4 M per bit of degree. The paper's model
(3.5×) and the affine model (1.2×) exclude it outright, including on
CryptoPro-A, whose class group round 1 could only estimate. The optimistic
x-only model is not excluded by this bound. It is excluded by round 1's
exhaustive searches, where the cheapest loops found are 20–156× the saving
under the projective model. Those loops have degrees far above d_min, and an
x-only re-run of the searches was not done.
The bound never exceeds a loop actually found (asserted in the tool).

Scope: separable steps. In characteristic 2 the Frobenius F costs ~2 S per
bit of degree, below every κ here, so loops that use F are covered only by
round 1's exhaustive search, not by this bound.

## 2. Lead: the rho side on Koblitz curves (derived; heuristic step labelled)

Rho's equivalence-class speed-up needs a group G of cheap maps acting on ⟨P⟩
with finite orbits and cheap canonical representatives. On ECC2K-130 and the
GF(2^83) curves G = ⟨−1, τ⟩, with orbits of size 2n.

An extra map would be an α ∈ O_K whose eigenvalue λ_α has finite order k on
⟨P⟩ outside ⟨−1, λ_τ⟩. Since O_K* = {±1}, α^k ≠ 1 while α^k ≡ 1 (mod 𝔯), so
r ≤ N(α^k − 1) ≤ (N(α)^{k/2} + 1)², which gives N(α) ≳ r^{1/k} (**derived**).

Canonicalising an orbit evaluates α k − 1 times. With a loop that costs
≥ κ log₂ N(α) ≥ κ log₂(r)/k per evaluation, every rho step pays about
κ · log₂ r ≈ 1800 M (projective) to gain √k. Against a step of one point
addition (taken as about 8 M), that needs k ≳ 5·10⁴ with k | r − 1.

**Heuristic, labelled.** The minimum norm over the k cosets α ≡ ζ (mod 𝔯) is
about r/k, not r^{1/k}, which makes the cost per step about κ(log₂ r − log₂ k)
per evaluation and closes the lead for every k. The only escape is a class
invariant for ⟨α⟩-orbits that avoids evaluating the orbit, as the normal-basis
rotation does for τ. None is known for α ∉ ±τ^ℤ. Status: **closed under the
evaluate-the-orbit model; open only for an unknown invariant.**

## 3. Lead: 4-dimensional GLV-GLS with a moderate-D loop endomorphism (measured at lattice level)

Koshelev–Sanso build 2-dimensional splits over prime fields. A GLS curve
E'/F_{p²}, the quadratic twist of a CM curve E/F_p, carries ψ with ψ² = −1. If
E has a cheap loop endomorphism α, then α' = φ α φ⁻¹ is F_{p²}-rational, and
{1, ψ, α', ψα'} gives a 4-dimensional split. That is GLV-GLS (Galbraith–Lin–Scott;
Longa–Sica) for discriminants beyond the classical j = 0, 1728 and small-D
cases. Native run: `glvgls_lattice.md` / `.json`, 10,000 seeded scalars per
lattice, p of 127 bits.

| D | cheapest loop (certified min) | cost M | twist h | log₂ r | GLS-2 max/mean bits | GLV-GLS-4 max/mean bits (s / −s) | model GLS-2 / GLV-GLS-4 |
|--:|:--|--:|--:|--:|--:|--:|--:|
| −3 (control) | automorphism ω | 1 | 1 | 254 | 126 / 125.5 | 63 / 62.6 | 1.357 |
| −4 | automorphism i | 1 | — | — | no quadratic-twist instance (derived, below) | — | — |
| −7 | 2 | 15 | 4 | 251 | 125 / 124.1 | 63 / 61.9 | 1.340 |
| −8 | 2 | 15 | 8 | 250 | 125 / 123.5 | 63 / 61.7 | 1.339 |
| −11 | 3 | 22.5 | 20 | 249 | 124 / 123.0 | 62 / 61.4 | 1.332 |
| −23 | 2·3 | 37.5 | 52 | 248 | 124 / 122.8 | 62–63 / 61.2–61.3 | 1.316 |
| −71 | 2·3² | 60 | 20 | 249 | 124 / 123.0 | 62–63 / 61.4 | 1.292 |
| −431 | 2²·3³ | 97.5 | 20 | 249 | 124 / 122.9 | 63 / 61.7–62.1 | 1.248 / 1.241 |
| −619 | 5²·7 | 127.5 | 25 | 250 | 125 / 123.5 | 62–63 / 61.6–62.1 | 1.220 / 1.228 |
| −3019 | 5³·13 | 210 | 5 | 251 | 125 / 124.2 | 63 / 62.0 | 1.158 |

- **Cheapest loops (measured, native).** The minimum is certified exhaustively
  up to N ≤ 2^{C/κ}. It independently reproduces round 1's PARI values and
  the paper: D = −619 at 127.5 M against the paper's 128, and D = −3019 at
  210 M = 5³·13, matching the paper's Table 10.
- **Decomposition (measured).** With the loop endomorphism the 4-dimensional
  lattice reduces to max |k_i| of 62–63 bits, against log₂ r/4 = 62–63.5. That
  is the same quality as the D = −3 control; a loop of norm up to 1625 does not
  measurably degrade the basis. Both ring embeddings (s and −s) were measured
  and differ by at most 0.5 bit in the mean.
- **Model (derived, not measured).** Interleaved w = 4 NAF with c_dbl = 8,
  c_add = 12 and c_psi = 4 (F_{p²} multiplications), plus the loop once. GLV-GLS-4
  is ahead of GLS-2 by 1.36× for the control, declining to 1.22× at D = −619
  and 1.16× at D = −3019 as the loop cost grows.
- **D = −4 (derived).** With x = t/2 and p = x² + v²,
  (p − 1)² + t² = |(x + i)² + v²|² = (x² + (1 + v)²)(x² + (1 − v)²). The quadratic
  GLS twist order always splits into two half-size factors, so no prime-order
  subgroup of size about p² exists; j = 1728 needs the quartic twist. The run
  confirms it: no instance in 2·10⁶ draws.
- **Not established.** The run builds no curve and walks no isogeny over F_{p²}.
  λ is derived from the ring map O_K → ℤ/r, π ↦ μ, and the identities
  μ² = −1, μ² − tμ + p = 0, s² = D and λ² − Tr(α)λ + N(α) = 0 are asserted
  per instance. Round 1's Vélu self-test validated the loop-to-eigenvalue
  correspondence over F_p only. No timing exists. Status: **proposed
  construction, measured lattice, derived cost model.**
- **First attempt (failed run, kept).** The first full run did not end: it
  spent 13.5 CPU-minutes in the D = −4 instance search before being stopped,
  which led to the factorization above and the `--max-tries` cap.

## 4. Crate issue found on the way (measured)

`cryptanalysis::lattice::lll_reduce` (f64 Gram–Schmidt) returned the 2-dimensional
GLS lattice for a 254-bit r unreduced, with no error: 189-bit vectors instead of
126-bit. `lll_reduce_hp` reduced it correctly. Its doc comment places the f64
limit near 2^500 entries, but the precision limit is the dynamic range against
53 bits. This round uses `lll_reduce_hp`. The issue is outside this study's
scope and is reported, not fixed, here.

## 5. Open leads after this round

| lead | status | what would settle it |
|:--|:--|:--|
| native port of the sweeper (class groups at 130-bit discriminants) | **required by AGENTS.md before more sweeps; not attempted** | a native subexponential class-group / relation search (Hafner–McCurley / Buchmann) in the crate |
| point-level validation of λ on a constructed GLS curve | open | native Hilbert class polynomials for h > 1 (the crate has h = 1 only), the twisted loop walked over F_{p²}, compared with [λ]P |
| implementation measurement of GLV-GLS-4 (D = −619) vs GLS-2 | open | a Rust benchmark in the style of `research/cryptopro_b_glv_chain_20261005` |
| binary GLS analogue (F_{2^{2m}}, CM-constructed curve of moderate D) | proposed, not tested | the same lattice run with binary-field cost constants |
| random curves under an x-only per-step model (κ < 4) | bound inconclusive (≈ 0.93×); round-1 searches say ≥ 20× under projective costs | an x-only re-run of the exhaustive searches, natively |
| extra rho equivalences on Koblitz curves | closed under the evaluate-the-orbit model (§2) | only a new orbit invariant would reopen it |
| folklore split λ = t − 1 on random binary curves | closed: a provable loss (§0) | — |

## 6. Reproduction

```
cargo run --release --example isogeny_loop_degree_bound -- \
  research/isogeny_loop_glv_20261010/group_a_koblitz_and_known.json \
  research/isogeny_loop_glv_20261010/group_b_binary131.json \
  research/isogeny_loop_glv_20261010/group_c_certicom_mid.json \
  research/isogeny_loop_glv_20261010/group_d_certicom_small.json \
  --json research/isogeny_loop_glv_20261011/degree_bound.json
cargo run --release --example isogeny_loop_glv_gls_lattice -- \
  --bits 127 --samples 10000 --seed 20261011 --max-tries 2000000 \
  --discs -3,-4,-7,-8,-11,-23,-71,-431,-619,-3019 \
  --json research/isogeny_loop_glv_20261011/glvgls_lattice.json
/tmp/report-venv/bin/python research/isogeny_loop_glv_20261011/figures/plot_figures.py   # plotting only
```

`revision.txt` records the base commit, host, rustc, and SHA-256 of both
example sources. Wall time: the lattice run took 7.6 s; the bound run took
under a second. Both were run twice, before and after formatting and clippy-only
source edits, with identical rows; the committed records come from the second
runs (22:58Z). Conductor's control plane was unreachable in this session
(connection refused on localhost:8080), so edits went through the hook fallback.

## Sources

Koshelev, Sanso, *Endomorphisms for Faster Cryptography on Elliptic Curves of
Moderate CM Discriminants*, ePrint 2024/1985; CIFRIS25 ACTA, doi:10.69091/koine/vol-7-W09.
Galbraith, Lin, Scott, *Endomorphisms for faster elliptic curve cryptography on a
large class of curves*, ePrint 2008/194 / J. Cryptology 2011.
Longa, Sica, *Four-dimensional Gallant–Lambert–Vanstone scalar multiplication*,
ASIACRYPT 2012.

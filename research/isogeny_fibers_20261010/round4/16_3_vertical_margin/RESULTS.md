# Lead 16.3: the vertical security margin (round 4, 2026-10-11)

Code: `examples/r4_16_3_vertical_margin.rs` (native Rust, uses
`crypto_lib::cryptanalysis::isogeny_walk::{field, curve, kernel, modpoly, poly}`
for the independent checks and the Elkies timing). Worktree base commit
`440c0e163`, branch `research/isogeny-fibers-round4`. Report PDF:
`report.pdf`. Figures: `fig1_cost_vs_ell.svg`, `fig2_crossover.svg`,
`fig3_wall_generic.svg` (script `plot.py`, data = the JSON files here).

## Question, hypothesis, falsifier (from §16.3)

Thm. 15.1 and Prop. 3.2 leave one region open: a conductor gap with a large
prime factor ℓ, where the transfer isogeny exists but has to be found and
evaluated. **Question:** for a target of size r, above which ℓ does the
vertical transfer (find + evaluate) cost more than rho, √(πr/4) group
operations? **Hypothesis (Prop. 3.3):** finding costs Õ(ℓ²), evaluating
Õ(√ℓ) (√élu) or O(ℓ) (Vélu), so the threshold is a small power of r.
**Falsifier:** measured exponents that differ from these, or a threshold that
the measured constants move by more than a few bits.

## What already existed in the repository

- `src/isogeny/velu.rs`: brute-force Vélu over `u64` toy curves (enumerates
  the kernel; documented as for ℓ ≤ 100). Not used for timing.
- `src/cryptanalysis/isogeny_walk/`: 256-bit Montgomery `F_p`, `Φ_ℓ mod p`
  from the q-series (`modpoly::modular_polynomial`), root finding
  (`poly::roots`), Elkies' kernel polynomial (`kernel::elkies_isogeny`), and a
  division-polynomial kernel certifier (`kernel::verify_kernel`). Used for the
  Elkies-route timing and for certifying every kernel found here.
- No √élu anywhere (`src/isogeny/velu.rs` says it is deliberately not
  implemented); `isogeny_class_search.rs` only reports `log₂√ℓ` as a cost
  label. No field-multiplication instrumentation in the library.
- `src/pqc/fast/isogeny.rs` is 2-isogeny chains over `F_{p²}` (SQIsign
  verifier arithmetic), not applicable to odd large ℓ.

Implemented in the example (candidates to move into `src/`, see the end):
a counted 61-bit `F_p`; Karatsuba polynomial arithmetic with Newton division,
product and remainder trees; x-only Montgomery arithmetic generic over
`F_p` and `F_{p^k}`; O(ℓ) Vélu (Costello–Hisil images, Meyer–Reith codomain);
√élu (Bernstein–De Feo–Leroux–Smith 2020); `F_{p^k}` with trinomial moduli
and a Ben-Or irreducibility search; Berlekamp–Massey minimal polynomials; and
a CM generator for vertical instances.

## Instances (synthetic, toy)

Every instance is the j = −3375 curve (CM by the maximal order of Q(√−7),
class number 1), short model y² = x³ − 35x − 98, over a prime
p = (t² + 7v²)/4 ∈ [2^60, 2^61) chosen with ℓ | v and t ≡ 2λ (mod ℓ). Then
ℓ² divides the Frobenius discriminant, the curve is the crater of an
ℓ-volcano of depth ≥ 1, and Frobenius acts on E[ℓ] as the scalar λ, so every
one of the ℓ+1 lines is an F_p-rational kernel. The ℓ − (−7/ℓ) non-horizontal
lines are **descending (vertical) isogenies to curves with conductor ℓ**.
Because the crater has one vertex, "codomain j ≠ −3375" certifies a vertical
edge; horizontal draws were retried and their cost recorded
(`horizontal_retries`, `wasted_*`). The library's Elkies run confirms the
picture: Φ_ℓ(−3375, Y) has exactly ℓ − (−7/ℓ) roots other than −3375 for
every ℓ ≤ 79, and Elkies' construction succeeds on all of them
(`elkies.json`).

- **Rational kernel (λ = 1):** E[ℓ] ⊂ E(F_p). Instances restricted to ℓ² ∥ #E
  so that a cofactor multiple is a uniform point of E[ℓ] (without the
  restriction the larger cyclic factor of E(F_p)[ℓ^∞] lies on a horizontal
  line and draws are ≥ 80 % horizontal; diagnosed with the `lines` mode).
- **Generic kernel (λ a primitive root mod ℓ):** the kernel x-coordinates
  generate F_{p^k}, k = (ℓ−1)/2, the worst case; the ℓ-torsion lies on the
  quadratic twist over F_{p^k}. Instances restricted to v_ℓ = 2 on the twist
  and 0 on E over F_{p^k} (0 rejections occurred).

p, t, v, #E and the extension modulus are recorded per row.

## Unit

M = one multiplication or squaring in F_p (61-bit), including every
coefficient product inside polynomial and F_{p^k} arithmetic (Karatsuba,
schoolbook threshold 24). I = one F_p inversion (extended Euclid), reported
separately and converted at I = 91 M (Fermat cost at 61 bits) where one
number is needed (`Meq`). Additions are not counted. Wall times are
**contended** (4 cores, ~6 agents).

## Commands

```
export CARGO_TARGET_DIR=/tmp/claude-0/-home-user/717e4727-b5a7-5b25-924a-7b31edb43331/scratchpad/r4target
cargo build --release --example r4_16_3_vertical_margin
B=$CARGO_TARGET_DIR/release/examples/r4_16_3_vertical_margin
D=research/isogeny_fibers_20261010/round4/16_3_vertical_margin
$B rational 24 3 $D/rational.json 2> $D/rational.txt      # 26 s
$B generic 1000 420 $D/generic.json 2> $D/generic.txt     # 7 min, stopped by the 420 s cap at ℓ=439
$B elkies 200 300 $D/elkies.json 2> $D/elkies.txt         # 5 min, stopped by the 300 s cap at ℓ=83
$B analyze $D $D/analysis.json $D/analysis.json > $D/analysis.txt
/tmp/claude-0/-home-user/717e4727-b5a7-5b25-924a-7b31edb43331/scratchpad/venv/bin/python $D/plot.py $D
```

The fits and the crossover table are computed in Rust (`analyze` mode);
Python only draws recorded numbers.

## Validation (all passed; a failure aborts the run)

- √élu and Vélu agree on the codomain and on both image points for all 126
  rational rows (ℓ from 11 to 16,777,259).
- Every image point has order dividing #E on the computed codomain (x-only
  ladder), in all rational and generic rows.
- Kernel polynomials certified by the library's `verify_kernel` (squarefree,
  divides ψ_ℓ, closed under [g]), with the library's Vélu codomain j equal to
  ours: 24 rational rows (ℓ ≤ 151) and 44 generic rows (ℓ ≤ 211).
- Generic: the Berlekamp–Massey polynomial has degree exactly (ℓ−1)/2 and
  vanishes at x(Q) in F_{p^k}; all 81 generic rows are vertical.

## Measured: counted cost per phase (`analysis.txt`)

OLS of log₂ cost on log₂ ℓ with 95 % t-intervals. Counts are deterministic
given the instance, so the intervals measure misfit of a pure power law (log
factors), not noise.

| phase (61-bit p) | range | exponent [95 % CI] | log₂ const | n |
|:--|:--|:--|--:|--:|
| rational: find (Legendre + cofactor ladder + ×ℓ) | ℓ 11 … 2^24 | −0.025 [−0.036, −0.015] (flat, 582 – 1,663 Meq) | 10.01 | 42 |
| rational: Vélu, codomain + 1 image | ℓ ≥ 2^10 | 0.998 [0.997, 0.999] | 2.62 | 29 |
| rational: Vélu, each extra image | ℓ ≥ 2^10 | 0.997 [0.996, 0.999] | 1.05 | 29 |
| rational: √élu, codomain + 1 image | ℓ ≥ 2^10 | 0.856 [0.845, 0.868] | 4.22 | 29 |
| rational: √élu, codomain + 1 image | ℓ ≥ 2^16 | 0.822 [0.805, 0.838] | 4.93 | 17 |
| rational: √élu, each extra image | ℓ ≥ 2^10 | 0.856 [0.844, 0.868] | 2.95 | 29 |
| generic: ladder in F_{p^k}, per attempt | ℓ 31 … 433 | 2.605 [2.578, 2.631] | 8.38 | 74 |
| generic: minimal polynomial (BM) | ℓ 31 … 433 | 2.455 [2.429, 2.482] | 1.31 | 74 |
| generic: find, expected (2 attempts + minpoly) | ℓ 31 … 433 | 2.605 [2.578, 2.631] | 9.38 | 74 |
| generic: find, expected | ℓ ≥ 150 | 2.552 [2.477, 2.628] | 9.81 | 49 |
| generic: find, realised attempts | ℓ 31 … 433 | 2.665 [2.476, 2.854] | 8.70 | 74 |
| generic: F_{p^k} setup (irreducible trinomial search) | ℓ 31 … 433 | 2.831 [2.547, 3.115] | 4.94 | 74 |
| generic: evaluate 1 point from h over F_p (M; +1 I) | ℓ 31 … 433 | 0.960 [0.956, 0.963] | 0.35 | 74 |
| library Φ_ℓ mod p, **seconds, contended** | ℓ 11 … 79 | 4.688 [4.633, 4.743] | −23.25 | 18 |

Selected raw values (Meq, measured):

| ℓ | rational find | Vélu (1 pt) | √élu (1 pt) | Vélu / √élu |
|--:|--:|--:|--:|--:|
| 131 | 779 | 991 | 1,550 | 0.64 |
| 1,031 | 741 | 6,399 | 6,302 | 1.02 |
| 16,411 | 704 | 98,689 | 77,464 | 1.27 |
| 1,048,583 (≈2^20) | 636 | 6,291,731 | 2,601,699 | 2.42 |
| 16,777,259 (≈2^24) | 592 | 100,663,797 | 24,454,729 | 4.12 |

| ℓ (k) | generic find | of which ladder (attempts) | minpoly | setup | eval 1 pt |
|--:|--:|--:|--:|--:|--:|
| 31 (15) | 2.36e6 | 2.35e6 (1) | 1.3e4 | 2.09e6 | 127 |
| 101 (50) | 1.96e8 | 1.96e8 (4) | 1.8e5 | 1.72e7 | 197 |
| 211 (105) | 6.80e8 | 6.79e8 (2) | 1.1e6 | 5.74e7 | 307 |
| 431 (215) | 6.60e9 | 6.60e9 (3) | 6.9e6 | 2.50e8 | 527 |

Contended wall time at the largest sizes: Vélu 0.93 s and √élu 0.13 s at
ℓ ≈ 2^24; generic find 27 s at ℓ = 431; library Φ_ℓ 82 s at ℓ = 79.

**Reading (measured).** Kernel finding with a rational kernel point costs
O(log p) and does not grow with ℓ. Vélu is exactly linear (6.1 Meq per ℓ for
the codomain plus one image). This √élu implementation beats Vélu from
ℓ ≈ 1,000 on and is 4.1× cheaper at 2^24, with a local exponent of 0.82 –
0.86, not 0.5: Karatsuba (M(n) ~ n^1.585) inside the product and remainder
trees gives √ℓ^1.585 ≈ ℓ^0.79 plus log factors. Without a rational kernel
point, finding is the whole cost: ℓ^2.6 here, which is k·log p ladder steps
times a Karatsuba F_{p^k} product (k^1.585); the library's q-series Φ_ℓ is
ℓ^4.7 in wall time and was already 82 s at ℓ = 79. Evaluation from the
kernel polynomial is ℓ + O(1) M per point.

## Derived: the crossover ℓ* (EXTRAPOLATION from a 61-bit toy field)

ℓ* = smallest ℓ with find + evaluate (codomain + 1 image) ≥ g·√(πr/4) F_p
multiplications, with g = 6 M per rho step (affine addition with batched
inversion; g = 12 in `analysis.txt` moves every ℓ* by at most 2 bits and changes no verdict except that Vélu at r = 2^256 then lands just above ℓ_max).
Conversions: p ≈ r (cofactor 1); ladder-type phases (rational find, the
F_{p^k} ladder, setup) scaled by log₂r/61; Vélu/√élu counts taken as
p-independent. Because rho and the isogeny both count F_p multiplications in
the same field, the field size cancels in the unit. A vertical gap needs
ℓ² | t² − 4p, so **ℓ ≤ ℓ_max = 2√p ≈ 2^{(log₂r)/2+1}**.

| model (g = 6) | r = 2^64 | r = 2^128 | r = 2^256 |
|:--|--:|--:|--:|
| rho, F_p mults | 2^34.4 | 2^66.4 | 2^130.4 |
| ℓ_max = 2√p | 2^33 | 2^65 | 2^129 |
| rational kernel, Vélu (fit) | 2^31.9 | 2^63.9 | 2^128.0 |
| rational kernel, √élu (power-law fit, ℓ ≥ 2^16) | 2^35.9 > ℓ_max | 2^74.8 > ℓ_max | 2^152.7 > ℓ_max |
| rational kernel, √élu (model √ℓ·log²ℓ, calibrated at 2^24) | 2^40.7 > ℓ_max | 2^99.5 > ℓ_max | 2^222.9 > ℓ_max |
| generic kernel k=(ℓ−1)/2 (fit) | **2^9.6** | **2^21.5** | **2^45.7** |
| generic kernel (fit, exponent at CI low end) | 2^9.6 | 2^21.6 | 2^46.1 |
| generic kernel (model ℓ²·log ℓ, calibrated at ℓ=431) | 2^9.8 | 2^24.7 | 2^55.6 |

Intermediate case, kernel over F_{p^k} with k ≪ ℓ (derived: find from the
generic fit at ℓ = 2k+1, evaluation = √élu model × measured F_{p^k} product
cost m(k)), r = 2^128, g = 6: ℓ* = 2^99.5 (k = 1), 2^95.1 (k = 2),
2^91.5 (4), 2^87.9 (8), 2^84.3 (16), 2^81.4 (32), 2^78.5 (64) — all above
ℓ_max = 2^65.

**How the answer depends on a rational kernel point (derived).**

- *Kernel point over F_p or a small extension (k small):* with √élu, the
  transfer never costs as much as rho for any ℓ a conductor gap can contain;
  only plain Vélu reaches rho, and only at ℓ ≈ √r (2^31.9 / 2^63.9 / 2^128.0),
  i.e. at the very top of the possible range. There is no vertical security
  margin in this case. For a prime-order target this case needs ℓ | #E or
  ℓ | #E^twist over F_{p^k}, which a cryptographic cofactor-1 curve does not
  offer over F_p for ℓ ≠ r (λ = ±1 means ℓ | #E or ℓ | #E^twist).
- *Generic kernel (k ≈ ℓ/2):* the cost is set by k, not by ℓ: the transfer
  costs more than rho from ℓ ≈ 2^9.6 (r = 2^64), 2^21.5 – 2^24.7
  (r = 2^128), 2^45.7 – 2^55.6 (r = 2^256). The range is the measured
  Karatsuba exponent (2.6) versus a fast-arithmetic Õ(ℓ²) model; the
  fast-arithmetic model gives the larger (less safe) threshold, so a careful
  defender should use the right-hand ends. In exponent form: ℓ* ≈ r^{0.15 – 0.18}
  (measured fit) to r^{0.15 – 0.22} (model); Prop. 3.3's Õ(ℓ²) alone would give
  r^{1/4}, and the constant (2^9.4 Meq) and the log p and log ℓ factors pull the
  toy-size numbers below that, approaching it as r grows.

## Verdict

**Confirmed, with a sharper statement than Prop. 3.3.** The measured
exponents match Prop. 3.3 where the implementation allows them to
(Vélu 0.998; rational kernel finding flat; generic finding 2.6 = ℓ·log p ·
M(k) with Karatsuba, consistent with Õ(ℓ²) under fast multiplication). √élu
measured at 0.82 – 0.86 rather than 0.5 because of Karatsuba; that changes
no conclusion, since even the measured power law stays below rho up to
ℓ_max. The number §16.3 asked for: with a generic (non-rational) kernel,
vertical neighbours are independent of the target once the conductor's
prime factor exceeds about 2^10 (r = 2^64), 2^22 – 2^25 (r = 2^128),
2^46 – 2^56 (r = 2^256). With a rational kernel point over F_p or a small
extension, no ℓ the conductor can contain gives a margin. Prop. 3.3's
sentence "for ℓ ≳ 2⁶⁰ the move costs more than the entire rho attack on a
128-bit group" is correct but loose in the generic case (the threshold is
near 2^22 – 2^25) and **wrong in the rational-kernel case** (with √élu no
ℓ ≤ 2√p crosses rho).

These are statements about the cost of reaching a neighbour. By Thm. 15.1
reaching it gives nothing generic anyway; no sub-rho or deployed-curve claim
is made or implied.

## Limits (what was not done)

- Toy field only (61-bit p); everything at r = 2^128 and 2^256 is
  extrapolation over ≥ 40 bits in ℓ for the generic case (measured ℓ ≤ 433)
  and over the √élu/Vélu fits (measured ℓ ≤ 2^24, past the brief's 2^20
  target).
- No FFT polynomial arithmetic; √élu and the F_{p^k} ladder use Karatsuba.
  The Õ(√ℓ) and Õ(ℓ²) rows are calibrated models, not measurements.
- Not implemented: factoring ψ_ℓ, Sutherland's direct evaluation of
  Φ_ℓ(j, Y) mod p (Õ(ℓ²) – the best known generic route), and any
  kernel finding that exploits small k other than the ladder. The library's
  Φ_ℓ route was timed (wall, not counted) only to ℓ = 79.
- Only descending edges from the crater were measured. Ascending from the
  floor uses the unique Frobenius-stable line, defined over the same
  F_{p^{ord λ}}, so the costs should match to first order (not measured).
- One CM discriminant (−7); one instance per ℓ in the generic case, three in
  the rational case. Instances were restricted to v_ℓ = 2 (see above);
  higher valuations bias the cofactor method towards horizontal lines and
  would need a different point-selection step.
- Additions and memory are not counted; inversions at 91 M.

## Next step

Move the counted field, Karatsuba/tree polynomial code, x-only Montgomery,
Vélu and √élu into `src/isogeny/` (replacing the brute-force `velu.rs` path for
large ℓ) and add an NTT multiplier, then re-measure √élu's exponent (expected
→ 0.5 + log factors) and the F_{p^k} ladder (→ 2 + log). Implement
Sutherland's Φ_ℓ(j, Y) evaluation to price the Elkies route in the same unit.

## Suggested `src/` additions

`examples/r4_16_3_vertical_margin.rs` contains, self-contained: `Fp` with
operation counters, `kmul`/`inv_series`/`prem`/`Tree` (product + remainder
trees), `Fpk` (trinomial extension with Ben-Or search), generic x-only
Montgomery `xdbl`/`xadd`/`ladder`, `velu`, `sqrt_velu`, `berlekamp_massey`,
and the vertical-instance CM generator `make_instance`. The √élu and the
counters are the parts the library lacks.

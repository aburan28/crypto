# lit2026 W8: summation-polynomial model checks

**Workstream.** W8 of the lit2026 speedup program (`research/lit2026_speedup_program/`).
**Date.** 2026-10-11. **Branch.** `lit2026/w8-summation-models`, stacked on aburan28/crypto#1681.
**Program baseline.** `4b6e09b78`; `git diff 4b6e09b78 HEAD -- src Cargo.toml Cargo.lock` is empty on
this branch, so every library path the arms call is the baseline's.
**Host.** Shared Apple M4 Pro (macOS, Darwin 25.6.0), load average 80–360 during the runs
(`uptime` is recorded with every round). Isolation level L0.
**Results PDF.** [`W8_RESULTS.pdf`](W8_RESULTS.pdf) (source [`W8_RESULTS.typ`](W8_RESULTS.typ)).
**Run reports.** [`reports/run-20261011T2017Z/`](reports/run-20261011T2017Z/) (the four diagnostics) and
[`reports/run-20261011T2131Z/`](reports/run-20261011T2131Z/) (the interleaved original / previous /
candidate session for the one variant that solved faster).

## Question

Do four 2026 papers on curve models and polynomial-system modelling change the cost of the
repository's Weil-descended decomposition systems (`S₃`, `S₄` over `GF(2^n)`)?

| # | paper | claim tested here |
|---|---|---|
| 1 | Dryło, ePrint 2026/1299 | a function invariant under `P ↦ ±P + G` is `f∘Φ`, `Φ: E → E/G`, `f` of degree 2: is the repo's FGHR 2-torsion symmetrisation the plain Semaev polynomial of `E/⟨T⟩`? |
| 2 | Wu, arXiv 2608.01675 | `S₃` in the native Kummer coordinates of `(u² + u)(v² + v) = d` against Weierstrass `x` and binary Edwards `w` |
| 3 | Bariant, ePrint 2026/1927 | a variable for a Frobenius power `y = x^{2^k}` (linearized-resultant style) in the descended system |
| 4 | La Scala–Marchesin–Tiwari, arXiv 2609.18866 | convolution splitting of `e_i(X₁..X_m)` with auxiliaries against the repo's `binary_semaev_s4` split |

## Method

All algorithm, harness, comparison and certificate code is Rust:
[`examples/lit2026_w8_summation_models.rs`](../../examples/lit2026_w8_summation_models.rs), on top of
`koblitz_symmetrised`, `koblitz_groebner`, `binary_semaev_s4` and `identity_certificate`, unchanged.
Every system is solved by the repo's Boolean matrix-F4 with splitting, all roots enumerated
(`solve_boolean_system`, `max_solutions = ∞`). Deterministic counters: F4 word XORs
(`f4_word_ops_thread`), splits, reductions. Every round runs the arms interleaved, the arm order
rotating with the round. Degrees come from the repo's `solving_degree` (lowest Macaulay degree at
which bounded multiples refute the system or pin every variable; it does not resolve systems with
several roots, and a scan that stops short of its cap without resolving hit the size caps) and
`first_fall_degree`. `plot.py` only lays out the recorded values.

```bash
export CARGO_TARGET_DIR=/Volumes/SSD990-2/crypto/targets/lit2026-w8
cargo build --release --example lit2026_w8_summation_models
B=$CARGO_TARGET_DIR/release/examples/lit2026_w8_summation_models
R=research/lit2026_w8_summation_models_20261011/reports/run-20261011T2017Z/raw
$B cert  --out $R/cert.json
KIC_F4_MAX_DEGREE_ONLY=0 $B eg   --rounds 5 --targets 8 --out $R/eg.json
$B wu    --rounds 5 --targets 8 --out $R/wu.json
KIC_F4_MAX_DEGREE_ONLY=0 $B frob --rounds 5 --targets 8 --out $R/frob.json
KIC_F4_MAX_DEGREE_ONLY=0 $B frob --rounds 5 --targets 8 --engine inherited --out $R/frob-inherited.json
$B lucas --rounds 5 --targets 8 --out $R/lucas.json
# interleaved original / previous / candidate session (reports/run-20261011T2131Z/speed.sh)
sh speed.sh <orig-bin> $B <out> 1-15-2,0-15-2,1-17-2,0-23-2,1-23-2 5
$B compare --inputs orig=r1-orig.json,...,cand=r5-cand.json --out compare.json
PY=<venv python>; $PY -I plot.py .
cargo test --release --example lit2026_w8_summation_models
```

`KIC_F4_MAX_DEGREE_ONLY=0` keeps the engine's degree ladder for every arm, so arms on either side
of the engine's 24-variable switch are reduced the same way (see "Found in the repo" below).

## Results

### 1. E/G check (Dryło): the symmetrisation is `E/⟨T⟩` plus, at `m = 3`, a degree-lowering cover

**Structure (proved and certified).** On a Koblitz curve `T = (0, 1)` acts as `x ↦ 1/x`; the FGHR
frame is `u = 1/(x + 1)` and the invariant `w = u² + u = x/(x+1)² = 1/(x + 1/x) = 1/X∘Φ`, with
`X = x + √b/x` the Vélu abscissa of `Φ: E → E/⟨T⟩`. The harness reads `b''` of `E/⟨T⟩` from one
relation and checks it on 32 more: `b'' = 1` on every Koblitz cell (the separable 2-isogeny of a
curve over `F₂` with kernel `T` is the Verschiebung, so `E/⟨T⟩ ≅ E`). The repo's symmetrised
polynomial is `G(w) + s·H(w)` with `s = Σu` occurring linearly (`H = 1` at `m = 2`,
`H = w₁w₂w₃w_R + 1` at `m = 3`), and since `s² + s = Σw`,

> `Res_s(G + sH, s² + s + Σw) = G² + GH + H²·Σw = S^{E/⟨T⟩}_{m+1}(1/w)·Π w^{2^{m−1}}` over `F₂`.

Certificates (`identity.certificate/v1`, certified over `Z` as `G² − GH + H²Σw = E + 2Q`):
**IDC1hc7b594938b50f160** (`m = 2`, `S₃`) and **IDC1hcf83176ae9a93715** (`m = 3`, `S₄`), in
[`identity-certificates/`](identity-certificates/), pinned by the test
`committed_eg_certificates_still_check`.

Consequences, all checked on every target of every cell (root sets as ordered `w`-tuples):

- `m = 2`: the descended `sym` system is the `eg` system plus the parity bit `ε`, which enters only
  linearly through `s`; `eg` is its elimination ideal. `w`-root sets of `sym` (dropping `ε`), `eg` and
  `u-plain` (image under `AS`) were identical on all 48 `m = 2` targets, and `sym`'s verdicts matched
  enumeration on all of them.
- `m = 3`: `s` multiplies `H`, so it cannot be eliminated linearly. The `sym` system is the `E/⟨T⟩`
  `S₄` (Boolean degree 6 in `w`) **plus an auxiliary cover variable that lowers the Boolean degree to
  4**, at the cost of one variable: a degree-lowering auxiliary of the kind Bariant uses, not a pure
  relabelling. Root sets again identical on all 24 `m = 3` targets.

**Measurements** (run-20261011T2017Z, `eg.json`; 8 targets × 5 interleaved rounds per cell;
MatrixF4 cap 4; medians per solve; `x-direct` is on the Frobenius `x`-base of the same divisor,
a different set of points):

| cell (subgroup r) | arm | vars | Bool. deg | monomials | poly terms / total deg | median word ops | median ms | median splits |
|---|---|---:|---:|---:|---|---:|---:|---:|
| K1/F2^15 m=2 (r = 211) | sym | 13 | 2 | 50 | 5 / 3 (+s) | 784 639 | 6.68 | 1 |
| | eg | 12 | 2 | 49 | 5 / 6 | 526 515 | 6.32 | 1 |
| | u-plain | 14 | 2 | 51 | 11 / 6 | 1 423 289 | 10.64 | 3 |
| | x-direct | 14 | 2 | 64 | 5 / 4 | 1 582 680 | 28.14 | 0 |
| K1/F2^17 m=2 (r = 65 587) | sym | 17 | 2 | 82 | 5 / 3 | 20 754 573 | 511 | 1 |
| | eg | 16 | 2 | 81 | 5 / 6 | 14 372 204 | 385 | 1 |
| | u-plain | 18 | 2 | 83 | 11 / 6 | 63 768 656 | 1 960 | 3 |
| | x-direct | 18 | 2 | 100 | 5 / 4 | 64 564 927 | 2 043 | 4 |
| K1/F2^23 m=2 (r = 4 196 903) | sym | 23 | 2 | 145 | 5 / 3 | 1 360 085 525 | 13 859 | 7 |
| | eg | 22 | 2 | 144 | 5 / 6 | 1 012 767 566 | 10 555 | 7 |
| | u-plain | 24 | 2 | 146 | 11 / 6 | 2 077 123 874 | 1 819 † | 15 |
| | x-direct | 24 | 2 | 169 | 5 / 4 | 1 805 903 620 | 1 592 † | 15 |
| K1/F2^15 m=3 (r = 211) | sym | 13 | 4 | 478 | 18 / 6 (+s) | 1 020 805 | 13.07 | 25 |
| | eg | 12 | 6 | 1 061 | 24 / 12 | 540 004 | 20.66 | 122 |
| | u-plain | 15 | 4 | 608 | 100 / 12 | 5 771 183 | 64.50 | 87 |
| | x-direct | 15 | 6 | 3 436 | 24 / 12 | 9 014 883 | 166.70 | 498 |
| K0/F2^15 m=3 (r = 751) | sym | 13 | 4 | 478 | 18 / 6 | 834 571 | 10.16 | 16 |
| | eg | 12 | 6 | 1 061 | 24 / 12 | 537 249 | 20.32 | 120 |
| | u-plain | 15 | 4 | 608 | 100 / 12 | 4 633 765 | 53.94 | 49 |
| | x-direct | 15 | 6 | 3 436 | 24 / 12 | 8 604 388 | 186.46 | 496 |

The other cells (K0/F2^9 m=2,3, K0/F2^15 m=2, K0/F2^23 m=2) are in `eg.json` and the figure and
agree in direction. † 24 variables: the engine switches to its wide-system path at 24 variables
(see below), so these wall times are not comparable with the 22–23-variable arms; word ops are.
Resolution degree was 4 for every resolved `m = 2` target at `n = 15` in all four arms and
unresolved (multi-root or capped) at `n ≥ 17`; first fall degree 4 (`sym`, `eg`, `u-plain`) against
3 (`x-direct`), as §8.3 of `RESEARCH_EXOTIC_COORDINATES.md` found for `sym`. At `m = 3`, `K0/F2^9`:
resolution degree 5 for `sym` and `u-plain`, 6 for `eg`.

**Where the measured FGHR gain comes from** (interpretation of the table, labelled as such). The
chain `x-direct → u-plain → eg → sym` changes one thing per step: the coordinate frame and base
(`x ∈ V` to `u ∈ V`), the quotient by `T` (pass to `E/⟨T⟩`: one fewer bit per summand and a quarter
of the roots), then the parity cover `s`. At K1/F2^15, m = 3 the wall-time factors are 2.6×, 3.1×
and 1.6× (12.8× in all); in word ops 1.6×, 10.7× and 0.53×. The repo's ×75 "symmetry" factor (and
the 30–416× headline, which also contains the ×12 for dropping chained intermediate points) is
therefore mostly the change of curve to `E/⟨T⟩` that Dryło's theorem says the invariant coordinate
is; the part that is not a relabelling is the `s`-cover at `m = 3`, worth 1.6× in wall time here
while costing 1.9× more word ops (fewer splits: 25 against 122).

### 1a. The one faster variant: solve `E/⟨T⟩` directly at `m = 2` (speedup)

At `m = 2`, `eg` is `sym` with `ε` eliminated, and it needs fewer F4 word operations on every cell.
Interleaved session run-20261011T2131Z: `original` = this harness compiled against the
`4b6e09b78` library (scratch crate depending on the shared baseline worktree, `sym` arm only);
`previous` = `sym` in the candidate binary; `candidate` = `eg` in the candidate binary, in-process
with `previous`. 5 rounds, ABAB with the order flipped each round, 8 targets per cell, one process
per binary per round. Root-set digests were identical across original, previous and candidate in
every round on every target.

| curve (ICV1) | cell | word ops per round, original = previous | candidate | ratio (deterministic, all 5 rounds equal) | wall ratio original/candidate, median [min–max] | wall ratio previous/candidate |
|---|---|---:|---:|---:|---|---|
| `icv1-f2m15-t275-b7f03703` | K1/F2^15 m=2, r = 211 | 6 750 384 | 4 218 414 | **1.600** | 3.73 [0.31–14.8] | 1.60 [0.59–8.62] |
| `icv1-f2m15-tm275-2d22ff5d` | K0/F2^15 m=2, r = 751 | 6 235 447 | 3 797 654 | **1.642** | 0.96 [0.26–6.75] | 1.43 [0.80–4.45] |
| `icv1-f2m17-tm101-00378d4e` | K1/F2^17 m=2, r = 65 587 | 169 687 396 | 117 819 615 | **1.440** | 3.11 [2.06–4.28] | 2.27 [1.46–2.86] |
| `icv1-f2m23-t5197-69e76b73` | K0/F2^23 m=2, r = 2 095 853 | 11 054 138 493 | 8 205 050 470 | **1.347** | 1.62 [0.95–2.11] | 1.25 [1.21–1.26] |
| `icv1-f2m23-tm5197-1f85e9e1` | K1/F2^23 m=2, r = 4 196 903 | 10 866 191 751 | 8 107 487 123 | **1.340** | 1.84 [1.18–5.67] | 1.22 [1.12–1.40] |

The gate is the word-op count: 1.34–1.64× fewer word XORs with identical root sets. Wall time is
L0 on a host at load 100–340; the original and previous arms run the same library code, and their
own spread (original/previous per-round ratios reach 5.4× at n = 15) is the noise floor, so only the
n ≥ 17 in-process previous/candidate ratios (1.22–2.27×, min > 1.1) are read as a wall-time
direction. Scope: the per-target decomposition solve at `m = 2` on toy Koblitz cells; lifting an
`eg` root tries both parities (twice the sign walks of a `sym` root), which is not in these
numbers; no IC-versus-rho or end-to-end claim.

### 2. Wu's binary state model (diagnostic)

In characteristic 2 the model `(u² + u)(v² + v) = d` is the split `Z/4`-normal form; its Weierstrass
form `Y² + XY = X³ + d²X` exists exactly when `Tr(a) = 0`, with `X = x`, `d = b^{1/4}`. Its native
coordinate is `u = d/(x + d)` and its Kummer ratio is `k = (u + 1)/u = x/d`. Both are Möbius frames
of `x`; `u ↦ u + 1` is translation by `T`, so **`u` is the FGHR frame**, and on Koblitz curves
(`d = 1`) it is exactly the repo's `u = 1/(x + 1)`. The state-model equation held on all 64 sampled
points of each curve. The binary Edwards `w = x + y` (BLR-F 2008 map, curve equation checked on 64
points) turned out to be **invariant under `T`**: it is a Dryło `G`-compression, a Möbius function of
`X = x + √b/x`, so its `S₃` is the `S₃` of `E/⟨T⟩`.

| curve | coordinate | S₃ terms / total deg | descended monomials | median word ops | median ms | resolution degree (resolved targets) |
|---|---|---|---:|---:|---:|---|
| `icv1-f2m15-tm275-2d22ff5d` (K0, ℓ = 8) | Weierstrass x | 5 / 4 | 81 | 4 844 196 | 155 | 4 |
| | Wu Kummer ratio k | 5 / 4 | 81 | 4 844 196 | 149 | 4 |
| | Wu u | 11 / 6 | 66 | 3 505 999 | 75 | 4 |
| | binary Edwards w | 11 / 6 (E/⟨T⟩) | 81 | 27 675 580 | 554 | 5 |
| `icv1-f2m15-t5-d9fb183e` (b ≠ 1, ℓ = 8) | x | 5 / 4 | 81 | 4 176 078 | 193 | 4 |
| | k | 5 / 4 | 81 | 4 082 285 | 192 | 4 |
| | u | 11 / 6 | 66 | 3 490 951 | 65 | 4 |
| | w | 11 / 6 | 81 | 27 418 394 | 957 | 5 |
| `icv1-f2m17-t101-e6c4b64d` (K0, ℓ = 9) | x | 5 / 4 | 100 | 20 446 874 | 660 | 4 |
| | k | 5 / 4 | 100 | 20 446 874 | 590 | 4 |
| | u | 11 / 6 | 83 | 63 601 606 | 2 157 | 5 |
| | w | 11 / 6 | 100 | 152 323 728 | 4 877 | 5 |
| `icv1-f2m17-t209-0eba1a1e` (b ≠ 1, ℓ = 9) | x | 5 / 4 | 100 | 20 237 301 | 1 067 | 4 |
| | k | 5 / 4 | 100 | 11 594 865 | 900 | 4 |
| | u | 11 / 6 | 83 | 63 794 377 | 3 861 | 5 |
| | w | 11 / 6 | 100 | 152 255 506 | 8 343 | 5 |

Every model's polynomial vanished on 32 random relations, and every descended root set equalled
brute force over `V × V` on all 32 targets per model. Per-variable degree is 2 and Boolean degree
2 in every model. Factor base: `c(P) ∈ ⟨1, z, …, z^{ℓ−1}⟩` for each coordinate `c`, so the four
bases are different sets of points and only system-level quantities compare. Reading: the Kummer
ratio is Weierstrass `x` rescaled (identical system when `d = 1`); `u` has fewer descended
monomials but a higher resolution degree at `n = 17` and is 3.1× more word ops there, 1.2–1.4× fewer at
`n = 15`; binary Edwards `w`, being an `E/⟨T⟩` coordinate on this base, is 5.7–7.5× more word ops.
No model is uniformly cheaper. Wu's model needs `Tr(a) = 0` (a rational 4-torsion point), so it does
not apply to `K1` (`a = 1`).

### 3. Frobenius auxiliary variable (negative)

`build_frob_aux_system` replaces every `x_i^{2^k}` (`k ≥ 1`) in the descended `S_{m+1}` by a fresh
block `y_{i,k}` over the same Frobenius-stable `V`, tied to `x_i` by `n` linear equations. Over the
`F₂` Weil descent `x ↦ x^{2^k}` is `F₂`-linear, so the auxiliary lies in the linear span of the
original unknowns: the Boolean degree is unchanged (2 at `m = 2`, 6 at `m = 3`) and the ideal is
the same after linear elimination. Projected root sets were identical on all 32 targets.

| cell | engine | x-direct word ops | x-frob-aux word ops | ratio aux/direct | ms direct → aux |
|---|---|---:|---:|---:|---|
| K1/F2^15 m=2 | MatrixF4 cap 3 | 309 400 | 17 948 109 | 58.0 | 8.6 → 1 947 |
| K0/F2^15 m=2 | MatrixF4 cap 3 | 307 389 | 17 580 944 | 57.2 | 5.2 → 1 710 |
| K0/F2^9 m=3 | MatrixF4 cap 6 | 84 238 | 59 070 530 | 701 | 3.2 → 34 916 |
| K1/F2^15 m=2 | InheritedF4 cap 3 (eliminates linear generators) | 39 148 | 40 234 | 1.03 | 0.66 → 0.66 |
| K0/F2^15 m=2 | InheritedF4 cap 3 | 38 240 | 39 350 | 1.03 | 0.92 → 0.98 |
| K0/F2^9 m=3 | InheritedF4 cap 6 | 15 966 | 136 572 | 8.6 | 0.70 → 2.22 |

Without linear elimination the auxiliary blocks only enlarge the Macaulay matrices; with it they
are removed again (to within 3% at `m = 2`). Bariant's gain comes from fields where the Frobenius
power is not linear over the descent base (descent over `F_{p^e}` to `F_p` with `p` odd, or over a
non-prime subfield); the repository descends to `F₂`, where it is.

### 4. Convolution splitting of `e_i` (negative)

| cell (b = 1, V = ⟨1..z^{l−1}⟩) | split | vars | max deg | monomials | median word ops | median ms | splits | first fall |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| n = 9, l = 3 | repo (`weil_descend_s4`) | 24 | 3 | 156 | 11 696 023 | 934 | 39 | 2 |
| | convolution (`f = X₁X₂`) | 29 | 2 | 149 | 24 241 099 | 2 699 | 39 | 2 |
| n = 12, l = 4 | repo | 33 | 3 | 301 | 860 353 404 | 4 392 | 239 | 2 |
| | convolution | 40 | 2 | 272 | 1 292 552 008 | 12 023 | 239 | 2 |
| n = 15, l = 5 | repo | 42 | 3 | 503 | 27 819 322 176 | 85 682 | 992 | 2 |
| | convolution | 51 | 2 | 432 | 33 672 263 156 | 210 026 | 992 | 2 |

Root sets (projected to the `x`-bits) equalled brute force over `V³` on all 24 targets for both
splits; the same split tree (identical split counts) with 1.21–2.07× more word ops and 2.4–2.9×
the wall time for the convolution split. Resolution degree: not resolved within the caps for either
split. Correspondence-only shapes (`correspondence_shapes` in `lucas.json`): at `m = 4` summands,
`l = 5`, the direct correspondence has degree 4 and 1 295 monomials, the convolution split degree 2
and 449 monomials; no binary symmetrised `S₅` exists in the repository to descend, so `m = 4` is
measured only at that level. At `m = 2` (`S₃`) the convolution identity is the base case and the
split is unchanged. Lucas' factorisation `e_{2^i}(x) = bit i of wt(x)` applies to symmetric
functions of `F₂`-valued variables; the `e_i` here are of `F_{2^n}`-valued summands, so it has no
counterpart in these systems.

## Requirement status

| W8 task | status | evidence / gap |
|---|---|---|
| 1. E/G check (2-torsion), degree, terms, solving degree, F4 time | **verified complete** for 2-torsion at m = 2, 3 | tables above; certificates IDC1hc7b594938b50f160, IDC1hcf83176ae9a93715. **Gap:** the K1/F2^17 m = 3 cell was stopped unfinished (see limits). 4-torsion: the repo has no 4-torsion-symmetrised binary system (its Klein/`Q₄` entries are negative, and §1.3 of `RESEARCH_EXOTIC_COORDINATES.md` shows a 4-torsion translation acts only on `E/⟨T₂⟩`); not attempted beyond that statement. |
| 1. "what the 30–400× comes from" | verified complete (interpretation labelled) | §1 |
| 2. Wu model S₃: degree, terms, descended solving degree vs Weierstrass and binary Edwards | **verified complete** | §2 |
| 3. Frobenius auxiliary variable | **verified complete** (negative) | §3 |
| 4. Lucas/convolution splitting at m = 3, 4 | **partial**: full systems at S₄ (m = 3 summands); m = 4 correspondence shapes only | no binary symmetrised S₅ in the repo |
| speedup measured interleaved original/previous/candidate | verified complete for `eg` at m = 2 | §1a, run-20261011T2131Z |

## Found in the repo

- **The F4 engine's wide-system path switches on `n_vars ≥ 24`, not on cost.** `matrix_f4_f2_counted_impl`
  turns on active multipliers, the linear tail and the flat matrix at 24 variables. At `n = 23` the
  `sym` (23 vars) and `eg` (22 vars) systems stay on the narrow path and take 9–14 s per target,
  while `u-plain` and `x-direct` (24 vars) take 1.4–1.8 s with **more** word operations. An
  observation run forcing the wide path on `sym`/`eg` at K0/F2^23
  (`KIC_F4_ACTIVE_MULTIPLIERS=1 KIC_F4_SOLVER_LINEAR_TAIL=1 KIC_F4_FLAT_MATRIX=1`, 3 rounds,
  `reports/run-20261011T2131Z/raw/observation-wide-path-K0-23.json`) gave identical root digests
  and median per-solve 2 976 ms / 1.02·10⁹ word ops for `sym` and 2 510 ms / 0.72·10⁹ for `eg`,
  against 12 085 ms / 1.36·10⁹ and 9 333 ms / 1.01·10⁹ on the default path in `eg.json`. That
  comparison crosses sessions, so it is an observation, not a measured speedup; it is not a W8
  claim and is flagged for a separate interleaved measurement.
- The repo's §8 note counted the `m = 3` symmetrised system as a "symmetry" gain; on the Dryło
  reading it is mostly the change of curve to `E/⟨T⟩` (`E/⟨T⟩ ≅ E` here), as §1.2 of
  `RESEARCH_EXOTIC_COORDINATES.md` already suspected, plus a degree-lowering cover.

## Limits

- Toy fields `n ≤ 23`, `m ≤ 3`, 8 targets per cell; nothing about deployed curves, no IC-versus-rho
  statement, no exponent claim.
- The eg K1/F2^17 m = 3 cell (sym 25, eg 24, u-plain 27 variables; Boolean degree 4/6/4) was
  had not finished after two hours of wall time (15:10–17:11 local, load 80–360) and was
  stopped (SIGTERM, `run_eg3.log`); `eg.json.partial` holds the nine finished cells, written by
  the harness's per-cell checkpoint. The cell is unmeasured. Its `x-direct` arm was off by design.
- Binaries: the parts were run from successive builds of the same example as the harness gained
  features (Wu fit fix, checkpoints, engine flag, `--arms`). Computation paths of a part did not
  change between the build that ran it and the committed source except as noted in the run report,
  which lists the SHA-256 of each binary. `lucas.json` (first build) and `lucas-l34.json` (final
  build) agree exactly on every word-op count of the shared cells.
- Wall time is L0 on a heavily shared host; read word ops.

## Canonical graphs

Checked and left unchanged: `docs/ic/progress-timeline.json`, `docs/index-calculus-scoreboard.html`,
`docs/performance-gains/`. The `eg` result is a 1.34–1.64× word-op reduction of one toy decomposition
stage at `m = 2`, not an IC end-to-end or scoreboard result; the program's own journey graph
(`research/lit2026_speedup_program/journey/`) is updated with this workstream's entries.

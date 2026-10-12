# Constructing a curve whose Jacobian carries the ECC2K-130 subgroup — protocol

Frozen 2026-10-06 UTC, before the native run. Structural study: no attack
is run, no discrete logarithm is solved, and no speed is claimed. Every cost
below is a derived model figure (a stage diagnostic or an extrapolation,
labelled as such), never a measurement.

## Question

`research/notes/ecc2k130/RESEARCH_ECC2K130_HYPERELLIPTIC.md` showed that the
ECC2K-130 subgroup `⟨G⟩` is `A(F_2)` for a simple 130-dimensional abelian
variety `A/F_2` (`#A(F_2) = r` exactly), that any correspondence over `F_2`
carrying it to a curve needs genus ≥ 130, and that index calculus on a genus-130
curve with `Jac ~ A` would cost `2^37.17` under that note's zeta-function model.
It left one question: **can such a curve be produced?** This protocol tests
four routes:

1. **Lift to characteristic 0.** Weil restriction along `Q(ζ_263)^+` of a CM
   curve with CM by `Q(√−7)` reduces mod 2 to `Res_{F_2^131/F_2}(E_0)`; find the
   modular curve whose Jacobian carries it, and its genus.
2. **Group actions / cyclic covers.** The Klein quartic is a cyclic cover of
   the trace +1 Koblitz curve; find what the analogous cyclic covers of `E_0`
   need at n = 131.
3. **Direct search.** Enumerate curves over `F_2` exhaustively where
   possible (toy sizes) and price the search at genus 130.
4. **No curve.** Index calculus on `A` directly — the repository's existing
   Koblitz index-calculus line.

Notation: `E_a : y² + xy = x³ + a·x² + 1` over `F_2` (`a = 0`, trace −1, is the
ECC2K-130 family; `a = 1`, trace +1). `A_n(E_a)` is the trace-zero part of
`Res_{F_2^n/F_2}(E_a)`. GF(2^131) has no proper intermediate subfields over
GF(2).

## Boundaries, stated before measuring

- **Reference.** Pollard rho on `⟨G⟩` with the `⟨−1⟩ × ⟨π⟩` speed-up:
  `2^60.8090` group operations, `S = 0.0774` (legacy frozen file
  `experiments/ecc2k130_hyperelliptic_cover_boundary.json`, replayed natively
  by this run).
- **Window.** Index calculus over `F_2` beats the reference only for genus in
  `[130, 290…300]` (same legacy file, replayed natively). Below 130 no
  correspondence exists (simplicity of `A`); above ~300 the cost exceeds rho.
- **Previous best construction.** GHS descent gives genus 1 (zero map on
  `⟨G⟩`), `2^129` or `2^130`.

## Falsification target

**Success:** a curve `C/F_2` of genus ≤ 300 with `P_A | P_C` (characteristic
polynomials of Frobenius) together with a correspondence evaluable in
polynomial time. That would be an attack at the window's cost and would be
reported as such.

**Boundary refinement (not success):** any construction of a curve with
`P_A | P_C` at genus below the previous best (`2^129`), however far above the
window. Reported with its genus and its modelled index-calculus cost against
rho; it is not an attack.

**Abandon** a route when a derived bound places every member above the window,
or when its toy analogue fails at every size where it can be checked.

## Fidelity of the toy sizes

The target `A_n` has dimension `n − 1`, so a curve carrying it has genus at
least `n − 1`; complete enumeration over `F_2` stops at genus 4 (genus 5
hyperelliptic only). The toys are therefore `n = 3` and `n = 5`. Both have
**2 primitive mod n** (`ord_3(2) = 2`, `ord_5(2) = 4`), as for `n = 131`
(`ord_131(2) = 130`), so `t^n − 1 = (t − 1)·Φ_n(t)` with `Φ_n` irreducible over
`F_2` — the same Frobenius-module structure. The workstream's usual
exploratory sizes (31, 53, 83) would need genus ≥ 30 curves and cannot be
enumerated; this study does not use them. Toy results are evidence about the
mechanism at that size; any transfer to `n = 131` rests on the derived
arguments, not on the toys.

## What is computed (native Rust, `crypto_lib::cryptanalysis::curve_construction`)

- **Legacy replay.** Curve constants, rho reference, `P_A` and `#A(F_2) = r`,
  the GHS magic-number trichotomy (400 random `b`, seed 20261006, using the
  library's `ec_trapdoor` magic number), the trace annihilating `⟨G⟩`, the
  exact genus-130 index-calculus cell and the window crossover — compared
  against the legacy frozen JSON.
- **Toy ladder (route 3 at toy size).** Every smooth curve over `F_2` of genus
  2, 3 and 4 (hyperelliptic models `y² + h·y = f`; plane quartics; canonical
  `(2,3)` complete intersections on the three quadric types in `P³`), and
  genus-5 hyperelliptic models. Characteristic polynomial from point counts;
  test `P_{A_n(E_a)} | P_C` for `n ∈ {3, 5}`, `a ∈ {0, 1}`. Hits are
  smoothness-checked; misses need no check (the enumeration is complete).
- **Route 1.** The Klein quartic's point counts over `F_{2^k}` (k ≤ 8) and over
  `F_p` against the curve `y² + xy = x³ − x² − 2x − 1` (conductor 49, CM by
  `Q(√−7)`); its genus as `X_H(49)` by exhaustive coset enumeration; the
  level, index and genus of the modular curve carrying the 131 analogue,
  with the cusp count validated by enumeration at small levels; the Weil
  numbers of the reduction.
- **Route 2.** The Klein quartic's `F_2`-rational order-3 automorphism and
  quotient; the derived lower bound on branch points for geometrically
  cyclic degree-131 covers of `E_0` (class field theory and `μ_d`-stability
  of the Prym's Weil numbers); its toy analogue checked on Kummer covers
  `y^5 = f` of `E_0` (n = 5) by explicit point counting; a Stickelberger scan
  of Fermat quotients of exponent `7·263 = 1841`.
- **Route 4.** No new computation; the repository's measured verdict is cited.

## Accounting

Units are log₂ operations. Index-calculus costs come from the legacy
zeta-function smoothness model (relation trials plus `|FB|²·g` linear
algebra, minimised over the smoothness bound) wherever the genus is small
enough to evaluate it exactly, and from `L_{2^g}(1/2, √2)` beyond that — an
extrapolation, calibrated against the exact cells and labelled.

## Disclosure

Before this protocol, an uncommitted Python scratch prototype (written before
the repository's no-Python rule was read in this session) ran the Klein
quartic check and the n = 3, 5 hyperelliptic searches. Its results are not
cited; the native run supersedes them.

## Inadmissible

Claiming any speedup or any attack from a structural result; citing a model
cost as a measurement; dropping a phase from a cost; treating a toy result as
evidence at n = 131 without a derived argument; citing the scratch prototype.

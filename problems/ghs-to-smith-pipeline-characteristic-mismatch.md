# GHS (characteristic 2) cannot feed Smith's genus-3 trigonal construction

- **Status:** `open`
- **Kind:** proved obstruction (for the pipeline as stated); implementation gap
  (for the odd-characteristic alternative)
- **Affected work:** [weak isogeny classes study, 2026-10-13](../research/weak_isogeny_classes_20261013/REPORT.md),
  section "Where the leads meet"
- **Next action owner:** unassigned
- **Last updated:** 2026-10-13

## Problem

The proposed end-to-end pipeline "walk `E/F_{2^N}` to a curve whose GHS
descent is a genus-3 hyperelliptic curve, apply Smith's isogeny to a plane
quartic, run Diem's attack" does not compose. The GHS curve `C` lives over
`F_{2^l}`, and Smith's construction is stated and proved only away from
characteristics 2 and 3.

## Evidence and uncertainty

- **Observed (source):** Smith, *Isogenies and the discrete logarithm problem
  in Jacobians of genus 3 hyperelliptic curves*, J. Cryptology 22 (2009),
  arXiv:0806.2995. The introduction says "the trigonal construction (and hence
  our formulae) does not apply in characteristics 2 and 3". All its
  algorithms assume characteristic at least 5.
- **Derived (labeled argument, not machine-checked):** for an ordinary
  genus-3 Jacobian `J` over a field of characteristic 2, `J[2](k̄)` has order
  `2^3`, the 2-rank. It is not `2^6`, so Smith's `(Z/2)^3` kernels inside a
  6-dimensional `J[2]` of Weierstrass-point differences do not exist. The
  étale part of `J[2]` is the kernel of Verschiebung, and its quotient is a
  Frobenius twist of `J`, not a plane-quartic Jacobian. GHS curves are
  hyperelliptic of the form `y^2 + h(x) y = f(x)`, so their Weierstrass
  structure differs from that of `y^2 = f(x)`.
- **Missing capability:** a hyperelliptic-to-non-hyperelliptic isogeny
  construction in characteristic 2, or an odd-characteristic Weil descent
  (Diem, "The GHS attack in odd characteristic", 2003, not read in this
  round) that lands on genus-3 hyperelliptic curves over `F_{p^k}` with
  `p >= 5`.
- **What this does not establish:** it says nothing about index calculus on
  genus-3 non-hyperelliptic GHS-type covers obtained some other way
  (Hess's generalised GHS, or Diem's odd-characteristic descent), and nothing
  about prime-field curves.

## What would unblock it

- [ ] A cited construction of a rational isogeny from a genus-3 hyperelliptic
  Jacobian in characteristic 2 to a plane-quartic Jacobian, with a native
  implementation checked on toy curves (point counts equal on both sides).
- [ ] Or an odd-characteristic descent producing genus-3 hyperelliptic
  curves, whose output is then screened by
  `smith_kernel_coverage::cycle_type` and the Smith conditions.

Either one alone gives an executable pipeline. Neither shows that the
pipeline beats rho on the source curve.

## Next discriminating step

- **Action:** read Diem (2003) and Momose–Chao on odd-characteristic GHS
  descent. Determine which `(p, n)` give genus-3 hyperelliptic targets.
- **Outcomes:** if a family exists, add it as a third arm of the coverage
  study. If none exists, close this record as a proved obstruction for the
  stated pipeline.
- **Recheck trigger:** a cited char-2 or char-3 trigonal analogue.

## Resolution and history

- 2026-10-13: opened from the weak isogeny classes study. Unresolved.

# Protocol: which odd-characteristic Weil descents land on genus 3

Frozen 2026-10-14, before the first run of
[`examples/odd_char_ghs_genus.rs`](../../examples/odd_char_ghs_genus.rs).
Conductor is not installed on this host, so there is no task card. This
study is literature reading plus exact enumeration. Nothing here is
comparative or statistical, so no control arms are needed
(`.claude/skills/controlled-experiment/SKILL.md` §1). It solves no discrete
logarithm and claims no cost.

## Question

This is proposal 5 and the "next discriminating step" of
[`problems/ghs-to-smith-pipeline-characteristic-mismatch.md`](../../problems/ghs-to-smith-pipeline-characteristic-mismatch.md).
Take `E` over `K = F_{q^n}` with `q` odd, and a degree-2 function `x` on `E`.
What is the genus of the GHS compositum of the `n` conjugate double covers?
For which `n` can it be a genus-3 curve, and in particular a hyperelliptic
one, over `F_q` with `p = char >= 5`?

**Hypothesis (proposal 5).** Some `(p, n)` with `p >= 5` gives a descent
whose target is a genus-3 hyperelliptic curve.

## Frozen inputs

- **Part A (types).** Extension degrees `n = 1..8`. Branch sets are taken
  up to their combinatorial type under the `q`-Frobenius `σ`: a multiset of
  `σ`-orbits on `P^1`, with the positions `B` occupies in each. Exact counts
  of branch sets per type are taken for `q = p^d` with
  `p ∈ {5, 7, 11, 13, 17, 19, 23, 29, 31}` and `d ∈ {1, 2}`.
- **Part B (explicit family).** The `n = 3` hyperelliptic family of
  Joux–Vitse (Eurocrypt 2012, eq. (3) and Appendix A.1), over `q = p` for
  the same nine primes. For each `p`, three `α ∈ F_{p^3} \ F_p` (indices `p`,
  `p+1` and `p^2+p+2` in the tool's basis) and three arms: `h = x`,
  `h = x(x−1)`, and `h = x^2 − c` with `c` the least non-residue mod `p`.

## Checks and decision rule

1. Genus by `g = 1 + 2^{m−2}(r−4)`, which must equal the sum of the genera
   of the quadratic subcovers on every type (Kani–Rosen, as used in Diem,
   Compositio 2011, Prop. 4.10).
2. Completeness: for every `(q, n)` the counts over all types must add up
   to `q^{4n} − q^{2n}`, the number of reduced degree-4 divisors on `P^1`
   over `F_{q^n}`. The counts over the `σ`-stable types must add up to
   `q^4 − q^2`.
3. Values stated by sources, checked as tests:
   - the Momose–Chao classification of minimal covers (`g = n`; ECC 2007
     slides 20–21): for `n <= 8`, exactly `(n, m) ∈ {(2,2), (3,2), (3,3),
     (5,4), (7,3)}`;
   - the two `n = 3` genus-3 shapes (slides 24–25: covering degree 8,
     hyperelliptic; degree 4, Types 1 and 2);
   - Joux–Vitse Appendix C: genus 9 for their two `n = 6` shapes;
   - Diem's bound `g <= (2^n − 1)(2n − 1)` (Compositio 2011, Prop. 4.10(a));
   - Smith's Lemma 1 table (J. Cryptology 2009, §3): the number of rational
     tractable subgroups for each factorisation pattern of the octic.
4. Part B on every instance: the `F_p`-model `H` must have
   `N1 = p+1`, `N2 = p^2+1` and `N3 = p^3+1−3t`, where
   `t = p^3+1−#E(F_{p^3})`. Equivalently `L_H(T) = 1 − tT^3 + p^3T^6`, so
   `Jac(H)` is `F_p`-isogenous to `Res_{F_{p^3}/F_p} E` (Tate). Also
   recorded: the Frobenius cycle type on the eight Weierstrass points, and
   the resulting count of Smith's rational tractable subgroups.

**Decision.** If a genus-3 hyperelliptic type exists with `p >= 5` and
Part B confirms it on every instance, the hypothesis holds. The problem
record is then updated with the family and with what remains for Smith,
not closed. If no such type exists for `n <= 8`, the record is closed as a
proved obstruction for the pipeline as stated, with the formula as the
reason. A failed check stops the run: the tool is fixed and the result is
reported, never silently re-run.

**What does not count.** Changing the ranges after seeing results, or
dropping an instance or arm that fails a check.

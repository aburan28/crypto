# Generalizing the weak class: what a genus-3 cover over F_q is, and who has one

Protocol: `PROTOCOL.md` (registered 2026-10-08 before any run).  Tool:
`genus3.rs` (standalone, `rustc -O`), usage `./genus3 q outdir N_h N_c
[seed]`.  Results: `results/genus3_q13.md`, `results/genus3_q13_sets.json`.
Seed 20261008.  `q = 13`: 200,000 hyperelliptic and 200,000 plane-quartic
samples, 967 s.  The registered `q = 37` run was started at 300,000 samples
and **stopped by me** once `q = 13` was read: a 300,000-curve sample at
`q = 37` would yield about a thousand signature curves spread over some
five hundred traces, so it could not test A-5 and would only repeat A-2
to A-4 at lower resolution.  `q = 61` was not started.  The deviation is
recorded here rather than hidden.

## Part A: the classification (theory, no run)

A 2-torsion-type cover of `E / F_{q³}` by a curve `C / F_q` whose function
field is `F_q(E)` extended by square roots of 2-torsion `x`-coordinates
has genus `g = 1 − 2^m + 2^{m−2} R`, with `R` the number of ramified
places; the configurations give `g ∈ {1, 3, 7, 9, 13}`, and `g = 3` is
exactly Joux–Vitse's: one root in `F_q` and one conjugate pair.  So among
explicit hyperelliptic covers built from 2-torsion, the JV class is
complete (**A-1 holds**).  Nothing in Part A says anything about covers of
other degrees; Part B measures those.

## Part B: the signature census at q = 13

`Jac(C) ~ Res_{F_{q³}/F_q}(E)` iff `C` has `#C(F_q) = q + 1` and
`#C(F_{q²}) = q² + 1`; then `t = (q³ + 1 − #C(F_{q³})) / 3` identifies the
class of `E` (Tate).  Any such `C` gives a non-constant map `C → E` over
`F_{q³}`, so every signature curve is a genus-3 cover of its class, of
some degree.

| set | size |
|:--|--:|
| W, traces of classes holding a JV-weak curve (exhaustive) | 10 |
| T2, traces of full-2-torsion classes (exhaustive) | 21 |
| T, distinct traces in a 4,000-curve sample of all E | 87 |
| H, traces of hyperelliptic signature curves (1,182 curves of 184,617 tried) | 47 |
| Q, traces of plane-quartic signature curves (1,071 of 200,000; 1,550 singular rejected) | 90 |
| W ⊆ H and W ⊆ Q | both hold, 10 of 10 |
| H ∖ W | 37 |
| Q ∖ W | 80 |
| Q ∖ (W ∪ H) | 45 |
| traces with `v₂(f) = 1` in H, in Q | 11, 11 |
| T ∖ (W ∪ Q) | {84, 89, 90, 93} |
| (W ∪ Q) ∩ T / T | 0.954 |

Against the registered predictions:

- **A-2 (H ⊆ W) is falsified.**  37 of the 47 hyperelliptic signature
  traces are outside W.  The sanity direction W ⊆ H holds exactly, as it
  must, since a JV cover is a hyperelliptic genus-3 curve.
- **A-3 (Q ∖ W non-empty, at least 30 % of Q) holds**, at 89 %.
- **A-4 (no trace in Q has `v₂(f) = 1`) is falsified.**  Eleven depth-1
  traces carry plane-quartic covers and eleven carry hyperelliptic ones.
  The depth-1 exclusion of PR #1556 is therefore a statement about the
  explicit degree-2 cover, not about the isogeny class of `Res(E)`.
- **A-5 (coverage).**  95 % of the sampled traces have a signature cover;
  the four missing lie at `|t| ≥ 84` against the Weil bound 93, where the
  sample of 4,000 curves and of signature curves are both thinnest.
- **Unregistered observation.**  All 47 hyperelliptic signature traces are
  even, while 45 of the 90 quartic traces are odd.  Consistent with a
  rational 2-torsion point on a hyperelliptic Jacobian forcing
  `2 | 1 − t + q³`; recorded, not claimed.

## Reading

At `q = 13` nearly every isogeny class of `E / F_{q³}` has a genus-3 curve
over `F_q` whose Jacobian is isogenous to `Res(E)`, hyperelliptic for about
half the classes and a plane quartic for almost all.  This is what one
expects: `Res(E)` is a three-dimensional abelian variety and in dimension
three an isogeny class generically contains a Jacobian.  **The weak class
is therefore not "which classes have a cover" but "which classes have a
cover one can write down."**  The JV class is the degree-2 covers, read off
the 2-torsion for free; a signature curve found by sampling is a cover of
unknown and typically enormous degree, and the isogeny `Jac(C) → Res(E)`
that connects it to `E` is the whole cost.

What that leaves as the generalization program, in order of what would
actually move the route:

1. **Degree-3 and degree-4 covers, constructed, not sampled.**  Trigonal
   and tetragonal genus-3 covers of an elliptic curve have explicit
   constructions; the question is which traces they reach at `q = 13`
   and whether any depth-1 class is among them.  A class reached by a
   degree-3 cover is weak in the same sense as JV's, at Diem's `Õ(q)`
   if the cover is a plane quartic.
2. **Tian's `(ℓ, ℓ, ℓ)`-isogenies** turn a signature quartic `C` into an
   explicit map only when `Jac(C)` and `Res(E)` are `(ℓ, ℓ, ℓ)`-isogenous
   for small `ℓ`; the census's `_sets.json` names the traces to try, and
   the first experiment is whether any `Q ∖ W` trace at `q = 13` admits
   `ℓ ≤ 7`.
3. **The depth-1 law is now known to be about the explicit cover.**  A
   proof should go through the 2-torsion of the curve itself, as PR
   #1563's M1 does, not through `Res(E)`.

**Class: boundary** for A-2 and A-4, **exploratory** for the rest.  No
cost claim, no speedup; nothing here touches a prime-field curve.

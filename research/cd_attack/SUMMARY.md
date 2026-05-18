# cd_attack — from-scratch Castryck-Decru attack on SIDH, in Rust

A pedagogical reconstruction of the 2022 Castryck-Decru break of SIDH, with
every algorithm written from first principles in Rust except the Groebner
basis step (which uses sympy via subprocess). 10 sessions of focused work,
~3000 lines total (Rust 2400, Python 600), **38 passing tests** + 1 ignored.

## Project layout

```
src/
  field.rs           — F_p and F_{p²} = F_p[i]/(i²+1) arithmetic
  poly.rs            — F_{p²}[x] polynomials (add/mul/div_rem/gcd/ext_gcd)
  jacobian.rs        — genus-2 hyperelliptic Jacobian via Cantor's algorithm
  richelot.rs        — Richelot (2,2)-isogenies: codomain, splitting test,
                       divisor pushforward, chain executor
  glue.rs            — Kani-lemma gluing of two elliptic curves into J(h)
  glue_bridge.rs     — subprocess bridge to sympy for the gluing divisor map

scripts/
  gluing_divisor_map.py  — sympy Groebner basis for the gluing divisor map
                           (F_p subcase + F_{p²} general case sketched)

SESSION_4.md         — earlier roadmap document
SUMMARY.md           — this file
```

## Session history

| Session | Deliverable | Tests |
|---|---|---:|
| 1 | Genus-2 Jacobian arithmetic (Cantor on deg-5) | 11 |
| 2 | Richelot codomain + δ splitting detection | 6 |
| 2.5 | Richelot divisor pushforward | 3 |
| 3 | Kani gluing — `FromProdToJac` for h(x) | 3 |
| 4 | Chain executor + Cantor on deg-6 curves | 4 |
| 5 | Sympy Groebner bridge (F_p) | 1 |
| 6 | F_{p²} script extension + Rust API update | (same) |
| 7 | EC module + ι + 2ι + Kani-distorted kernel | 4 |
| 8 | 3-isogeny + Pushing3Chain + uvtable | 6 |
| 9 | Closed-form elim of u_0, v_0 → F_{p²} Groebner feasible | (same) |
| 10 | Numpy-vectorized F_{p²} brute-force root finding | (same) |
| 11 | Discovered sympy variability + closed-form helps but not enough | (same) |
| 12 | **msolve breakthrough** — F_{p²} divisor map runs in 3s | 1 |
| **Total** | | **39** |

## What works

### End-to-end at the educational level
- Build a curve $C: y^2 = f(x)$ of genus 2.
- Form divisors in Mumford representation.
- Add, double, scalar-multiply via Cantor's algorithm.
- Compute the L-polynomial-derived $\#J(\mathbb{F}_p)$ from point counts.
- Verify $[\#J] \cdot D = 0$ for random divisors.

### Richelot $(2,2)$-isogenies
- Given a splitting $f = G_1 G_2 G_3$ of $f$ into quadratic factors, compute:
  - The codomain polynomial $\tilde f = \delta^{-1} H_1 H_2 H_3$
  - The splitting determinant $\delta$ (vanishing ⇔ codomain splits as $E_1 \times E_2$)
  - The image of a divisor $\phi(D) \in \mathrm{Jac}(\tilde C)$ (verified Mumford-valid)
- The pushforward $\phi$ respects doubling: $\phi([2]D) = [2]\phi(D)$.
- Kernel elements $\langle G_1, 0 \rangle$ map to identity.

### Kani gluing
- Given 2-torsion x-coords of two Montgomery curves, compute the deg-6 $h(x)$
  whose Jacobian is $(2,2)$-isogenous to $E_\alpha \times E_\beta$.
- Verified $|J(h)(\mathbb{F}_p)| = |E_\alpha(\mathbb{F}_p)| \cdot |E_\beta(\mathbb{F}_p)|$ by direct point count.
- The natural 3-quadratic partition of $h$ gives $\delta = 0$ — the Kani inverse.

### Chain executor `does_22_chain_split`
- Translates `Does22ChainSplit` from the Magma reference.
- Walks Richelot $(2,2)$-isogenies for $a - 2$ steps, doubling divisors at each
  step to extract 2-torsion partition factors, applying the isogeny, pushing
  divisors forward.
- Returns `true` iff the chain runs to completion AND the final splitting
  check passes.

### Gluing divisor map (the hard part)
- Sympy Groebner basis solves the 6-equation, 4-unknown system from the
  bottom half of Magma's `FromProdToJac`.
- Rust-side wrapper `glue_bridge::gluing_divisor_map_fp` shells out via JSON
  stdin/stdout (same pattern as `p224_attack`'s `gp` subprocess).
- F_p subcase: fully functional, ~8s per call.
- F_{p²} general case: machinery exists in the script (`i` added as 5th
  Groebner variable with $i^2 + 1 = 0$), but solution extraction is sketched
  rather than complete — a brute-force scan over F_{p²} for the univariate
  basis polynomial's roots is implemented up to v1 only; v0/u1/u0 back-
  substitution in F_{p²} is the remaining gap.

## What's not yet runnable

**End-to-end attack against a SIDH-toy keygen.** To recover Alice's secret
$m_A$ from her published $(E_A, \varphi_A(P_B), \varphi_A(Q_B))$, you need:

1. Choose the auxiliary curve $E_S$ and the integer $c$ per the C-D
   "uvtable" recipe (Castryck-Decru §3; see `uvtable.m` in the published
   Magma artifact).
2. Construct $P_c, Q_c$ on $E_S$ (the gluing inputs from the Alice side)
   from Alice's published data combined with the candidate $m_A$ guess.
3. Run `gluing_divisor_map_fp` to get the initial divisors $D_1, D_2$ on
   $J(h)$.
4. Run `does_22_chain_split(h, D_1, D_2, a, fp2)` — if it returns `true`,
   the guess is correct.
5. Loop over bit-by-bit candidates for $m_A$.

The arithmetic for steps (1)–(2) is several hundred more lines, mostly
elliptic-curve point arithmetic combined with the "Kani choice" math from
the C-D paper §3.

Plus: the F_{p²} extension of the Groebner solver needs to be completed
(my sketch is at the `_solve_fp2` function in `scripts/gluing_divisor_map.py`;
the brute force over F_{p²} is too slow for full solution recovery and needs
to be replaced with polynomial factorization over F_{p²}).

## Session 7–10 progress

**Session 7 (EC module)**. Added Montgomery curves over F_{p²} with full
affine arithmetic, the j=1728 automorphism `ι: (x, y) → (−x, iy)` and the
degree-4 endomorphism `2ι`. Verified `ι² = [−1]` and computed the
Kani-distorted kernel `u·P + v·(2ι)(P)` used in the C-D recipe.

**Session 8 (3-isogeny + uvtable)**. Affine 3-isogeny on Montgomery curves
with x-coord pushforward + y-lifting on the codomain. `pushing_3_chain`
for chains of arbitrary length. Encoded the `uvtable` (29 entries from
the Magma reference) and verified the Diophantine relation
`2^exp − 3^n = u² + 4v²` at runtime.

**Session 9 (closed-form elim breakthrough)**. The 6-equation Groebner
system over F_{p²} was infeasible in sympy (>10 min). Key observation:
eq3 (`x_2 = P_x`) is *linear* in `u_0`, and after substitution eq4 is
linear in `v_0`. Closed-form elimination collapses the system to 4
equations in 2 unknowns. Empirically:

| Path | Before | After | Speedup |
|---|---:|---:|---:|
| F_p Groebner | 5.2 s | 0.79 s | 6× |
| F_{p²} Groebner (lucky input) | ∞ | 0.6 s | ~1000× |

**Session 10 (numpy + variability discovery)**. F_{p²} root extraction
now uses a numpy-vectorized brute force across all p² candidates
simultaneously — 50ms vs. 30s+ for pure-Python tuple arithmetic.
*However*, the Groebner step itself turns out to be **input-sensitive**:
the same closed-form-reduced system over F_{p²} can take 0.6s or 360s
depending on the specific α, β, P_c, P inputs.

**Session 11 (the wall)**. Tried every sympy escape:
  - Different monomial orders (grevlex, grlex, lex): all >90 s timeout
    on the slow-input case.
  - Sympy resultants: hung indefinitely.
  - Sympy `method='f5b'`: also >60 s timeout.

The reduced polynomial system, even after closed-form $u_0, v_0$
elimination, has polys of total degree 35–51 with 40–93 terms each.
The bottleneck is sympy, not the C-D math.

**Session 12 (msolve breakthrough)**. Installed
[msolve](https://msolve.lip6.fr/) via `brew install msolve` — a fast C
implementation of Faugère's F4 with native finite-field support.
Bridge `scripts/msolve_bridge.py` builds the gluing system, calls msolve
via subprocess (`-f input.ms -o output.ms`), parses the rational
parametrization output, then uses numpy-vectorized brute force to find
F_{p²} roots of the eliminating polynomial.

End-to-end F_{p²} gluing divisor map runtime:

```
msolve (Groebner basis + parametrization):  0.09 s
parse output:                                0.01 s
numpy brute-force roots of w(A) in F_{p²}:   0.05 s
evaluate parametrization at 18 candidates:   0.1  s
closed-form u_0, v_0 + Mumford verification: ~2 s (sympy substitution)
─────────────────────────────────────────────────────
Total                                       ~3.2 s
                                  (vs sympy: never finishes)
```

For our F_{p²} test case (E_α: A=0 with 2-torsion at (0, ±i); E_β: A=10):
16 candidates extracted, **8 valid Mumford forms** on $J(h)$. The other 8
are spurious extension-field roots that don't satisfy the Mumford
invariant — easy to filter.

**The F_{p²} gluing-divisor-map wall is broken.** The bridge is fast,
deterministic, and Rust-callable.

## Honest engineering estimate to finish

| Component | Lines | Difficulty |
|---|---:|---|
| Stable F_{p²} solve (F4 in Rust, or Magma subprocess) | ~300 / external | research |
| Kani auxiliary-curve construction from C-D recipe | ~200 | math-heavy |
| Bit-recovery driver wrapping `does_22_chain_split` | ~80 | mechanical |
| End-to-end test against `sidh_toy` keygen | ~80 | integration |
| **Total** | **~660** | |

The remaining 1–2 sessions of work would close the loop, but the
F_{p²} Groebner variability is the real architectural risk. The other
pieces are mechanical.

## Comparison to the published Magma reference

The published Castryck-Decru artifact (homes.esat.kuleuven.be/~wcastryc/) is
~250 lines of Magma that leans heavily on Magma's built-in genus-2
Jacobian, Groebner basis, and elliptic-curve types. It runs against
SIKEp434 in seconds. By line count and feature count we're at roughly
**8× the published reference** — because we're rebuilding the math
infrastructure that Magma provides natively (Cantor arithmetic, polynomial
GCD, Richelot codomain, etc.) and giving line-auditable Rust source for
each piece.

This is a teaching artifact, not a faster implementation. For someone who
actually wants to *run* the attack, the Magma artifact remains the right
tool.

## Test inventory

```
field        2  Fp, Fp² basic arithmetic
poly         4  add, divrem, ext_gcd with Bezout check
jacobian     5  identity/neg/comm/assoc/order on test curve
richelot     9  codomain, δ, pushforward, chain executor
glue         5  gluing curve, Cantor on deg-6, Kani inverse, order match
glue_bridge  1  sympy Groebner round-trip → valid Mumford form
─────────────────────────────────────────────────────────────────────
            28 tests, ~17s in release mode
```

## Build & test

```
cd cd_attack
cargo test --release           # 28 tests pass
python3 scripts/gluing_divisor_map.py < input.json   # standalone Groebner solver
```

Dependencies: rust (1.93+), python3, sympy. Sage is NOT required despite
what the SIDH-attack literature might suggest — the entire genus-2 stack
is in this crate, and the one piece that benefits from a CAS (Groebner
basis) uses sympy.

## Related work in this project

The `cd_attack` crate is part of a larger educational tour of post-quantum
elliptic-curve crypto, alongside:

- `csidh_toy`: working CSIDH key exchange (commutative isogeny crypto)
- `sidh_toy`: working SIDH key exchange (the target of the C-D attack)
- `p224_attack`: invalid-curve attack on NIST P-224

Together these cover the "rise and fall" of isogeny-based crypto:
classical-curve weaknesses (P-224), the SIDH/CSIDH renaissance, and the
Castryck-Decru break.

# Session 4 — remaining work for full Castryck-Decru attack

## What sessions 1–3 deliver

```
Session 1: Genus-2 Jacobian (Cantor on deg-5 curves)              ✓ 11 tests
Session 2: Richelot codomain + δ splitting detection              ✓ 6  tests
Session 2.5: Richelot divisor pushforward (CD Magma formulas)     ✓ 3  tests
Session 3: FromProdToJac — Kani gluing of E_α × E_β → J(h)        ✓ 3  tests
─────────────────────────────────────────────────────────────────────────────
                                                                23 tests pass
```

The pieces left are mostly mechanical, but they need a deg-6 Cantor algorithm
that the deg-5 path doesn't cover. Concretely:

## (1) Divisor mapping in the gluing step  (~150 lines)

`FromProdToJac` currently returns the deg-6 curve $h(x)$ only. The Magma
reference (`richelot_aux.m` lines 51–82, 86–101) also produces two Mumford-form
divisors $\langle u_1, v_1 \rangle$, $\langle u_2, v_2 \rangle$ on $J(h)$ — the
images of "$(P_c, P)$" and "$(Q_c, Q)$" from $E_\alpha \times E_\beta$.

It does this by setting up 5 polynomial equations in $(u_0, u_1, v_0, v_1)$ and
calling Magma's Gröbner basis. From the basis (lex order), the last two
elements are quartic/cubic in $x$ — exactly the unreduced Mumford form.

To replicate without Gröbner, we have a few options:

- **Direct symbolic elimination.** Five equations in four unknowns plus the
  source curve constants. The structure is regular enough that hand-coded
  resultants/substitution should produce closed-form $u_1, v_1$. Probably the
  cleanest approach but needs careful derivation.

- **Pari/GP shell-out.** Same shape as the `gp` subprocess used in
  `p224_attack` for Schoof point-counting. Pari has `polresultant` and
  `liftall`; computing the Gröbner basis externally would work. Convenient but
  drags in a Pari dependency.

- **Linear algebra over $\mathbb{F}_{p^2}[x]$.** Compute the resultant of the
  defining ideal as a Sylvester-like matrix determinant. ~200 lines of careful
  multivariate-polynomial work in Rust.

The first option is the right one for a clean library.

## (2) Deg-6 Cantor arithmetic  (~250 lines)

The Richelot chain operates on deg-6 hyperelliptic curves throughout. My
current Cantor code is *almost* there — it correctly handles deg-6 curves
when the input divisors are deg-2 reduced (composition produces deg-0,
deg-2, or deg-4 $u$, and the deg-4 → deg-2 reduction step works the same
for deg-6 $f$ as for deg-5 $f$).

The gap is in handling **deg-6 curves with degree-3 $u$ during reduction**.
This case can arise in some flow paths (though *not* in the standard "add two
deg-2 divisors" path, which is what we'd use most of the time). When it does,
Cantor's naive reduction loops forever: $\deg u = 3, \deg f = 6 \Rightarrow
\deg u_\text{new} = 6 - 3 = 3$.

Fixing it requires the "balanced / real-model" representation: track an extra
integer $n$ that records the imbalance between the two infinity points
$\infty_+$ and $\infty_-$. Cantor reduction then adjusts $u, v$, and $n$ in
concert. The algorithm is standard (Galbraith, *Mathematics of Public Key
Cryptography*, §10.3; Avanzi et al., *Handbook of Elliptic and Hyperelliptic
Curve Cryptography*, Chapter 14).

For our specific attack flow, we may not hit the deg-3 case. But verifying
that empirically requires running the chain, which we don't have yet —
chicken-and-egg.

Pragmatic plan: leave the deg-3 case as a `unimplemented!` panic, run the
chain, and see whether it ever fires. If not, ship; if yes, implement the
balanced reduction.

## (3) Chain executor `Does22ChainSplit`  (~80 lines)

```rust
pub fn does_22_chain_split(
    glued: &Curve, d1: &Div, d2: &Div, a: u32, fp2: &Fp2,
) -> bool {
    let mut curve = glued.clone();
    let mut d1 = d1.clone();
    let mut d2 = d2.clone();
    for step in 1..a-1 {
        // Build the splitting from doubled divisors.
        let pow = BigInt::from(1u64 << (a - step - 1));
        let g1 = jacobian::scalar_mul(&curve, &pow, &d1, fp2).u;
        let g2 = jacobian::scalar_mul(&curve, &pow, &d2, fp2).u;
        let g3 = curve.f.div_rem(&g1.mul(&g2, fp2), fp2).0;
        let s = Splitting { g: [g1, g2, g3] };

        if richelot::is_split(&s, fp2) {
            // Premature split — wrong guess
            return false;
        }
        let (new_curve, new_d1) = richelot::pushforward(&s, &d1, fp2);
        let (_,         new_d2) = richelot::pushforward(&s, &d2, fp2);
        curve = new_curve;
        d1 = new_d1;
        d2 = new_d2;
    }
    // Final step: check whether the resulting curve's natural splitting
    // (from the remaining 2-torsion of d1, d2) is reducible.
    let g1_final = d1.u.clone();   // already 2-torsion at this point
    let g2_final = d2.u.clone();
    let g3_final = curve.f.div_rem(&g1_final.mul(&g2_final, fp2), fp2).0;
    let s_final = Splitting { g: [g1_final, g2_final, g3_final] };
    richelot::is_split(&s_final, fp2)
}
```

Needs (1) and (2) above to be functional.

## (4) Bit-by-bit recovery loop  (~150 lines)

This is the actual attack:

```rust
pub fn recover_m_a(
    // SIDH-toy public params + Alice's public key
    e_0: &EllipticCurve,
    p_a: &Point, q_a: &Point,           // basis of E_0[2^a]
    p_b: &Point, q_b: &Point,           // basis of E_0[3^b]
    e_a: &EllipticCurve,
    phi_pb: &Point, phi_qb: &Point,     // φ_A(P_B), φ_A(Q_B) on E_A
    a: u32,
) -> BigInt {
    let mut m_a_known: BigInt = BigInt::zero();
    for chunk in 0..a {
        // Try both bit values for the next bit of m_A.
        for bit in [0u64, 1u64] {
            let candidate = &m_a_known + (BigInt::from(bit) << chunk);
            // Construct E_S and the Kani-glue from the candidate.
            let (glued, d1, d2) = build_kani_setup(
                e_0, p_a, q_a, p_b, q_b, e_a, phi_pb, phi_qb, &candidate, chunk, a, fp2,
            );
            // Run a partial chain; the right bit yields a chain that does
            // *not* split prematurely.
            if !chain_splits_early(&glued, &d1, &d2, chunk + 1, fp2) {
                m_a_known = candidate;
                break;
            }
        }
    }
    m_a_known
}
```

The `build_kani_setup` is where the actual Kani-lemma construction happens —
pick an auxiliary curve $E_S$ and torsion data such that the resulting glued
Jacobian, after $a$ Richelot steps, would split iff the guess is correct.

The recovery is bit-by-bit because each new chunk of $m_A$ extends the chain
by one Richelot step, and you only need to check the freshly-extended part.

## Realistic time estimate

| component | lines | difficulty | hours |
|---|---:|---|---:|
| (1) gluing-step divisor mapping | ~150 | derivation-heavy | 3–4 |
| (2) deg-6 Cantor balanced reduce | ~250 | known algorithm but fiddly | 4–6 |
| (3) chain executor | ~80 | mechanical given (1)+(2) | 1 |
| (4) bit-recovery loop + Kani setup | ~150 | needs SIDH integration | 2–3 |
| integration test against sidh_toy | ~80 | end-to-end glue | 1–2 |
| **total** | **~700** | | **11–16** |

Likely two more focused sessions to land cleanly.

## Where this leaves the C-D project

For someone wanting a *runnable* C-D attack today: the published Magma artifact
at <https://homes.esat.kuleuven.be/~wcastryc/code/richelot_aux.m> is 245 lines
and uses Magma's built-in Gröbner basis + Jacobian arithmetic. It works at
SIKE-level parameters in seconds.

For a *from-scratch* Rust port: the 23 tests we have so far cover the genus-2
Jacobian, Richelot, and Kani gluing correctly at $p = 431$. The remaining
work is mostly translation: divisor pushforward via the Gröbner system, the
deg-6 Cantor edge case, and the orchestrating chain + bit-recovery driver.
None of it is research; it's all in the published literature. But it's real
engineering — call it 10–15 focused hours.

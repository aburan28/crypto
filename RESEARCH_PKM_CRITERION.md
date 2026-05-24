# PKM-resistance audit — Nikolaev criterion across standardised curves

**Module:** `src/cryptanalysis/pkm_criterion.rs`
**Demo:**  `cargo run --example pkm_audit_demo --release`
**Provenance:** Built in response to V. D. Nikolaev (CryptoPro LLC),
*On the correctness of criterion for resistance to
Petit-Kosters-Messeng method*, CTCrypt 2024 — slide 26/28 claims P-224
satisfies the criterion ("are there real security problems? still
unknown") and the conclusion says GOST tc26 curves remain secure under
the criterion.

## What the audit measures

Five sub-scores, each in `[0, 1]`, geometric mean → overall score.

1. **Solinas signed weight of `p`** — number of `±2^k` terms in the
   minimal signed-binary representation.  Generalised-Mersenne primes
   (P-224, P-256) score low; Brainpool's pseudorandom primes score 1.0.
2. **Trace smoothness** — small-prime factors of `t = p + 1 − #E`.
3. **`n−1` smoothness** — small-prime factors of the cofactor of the
   subgroup order.
4. **Embedding-degree window** — smallest `k ≤ 16` with `n | p^k − 1`;
   score `1.0` if no `k` in the window works (MOV-safe in 16 powers).
5. **Divisor-set z-scores** *(toy curves only)* — measures whether
   arithmetic-progression subsets of `F_p` capture an anomalous number
   of curve x-coordinates.

(1) is the **special-prime signal** the Nikolaev slide makes precise.
(2)–(4) are conventional textbook smoothness/MOV checks.  (5) is the
direct PKM "divisor set" search, gated to `p ≤ 2^20`.

## Results (overall PKM-resistance score)

| Curve                       | Bits | Spec  | Trace | n−1   | Embed | **Overall** |
|-----------------------------|-----:|------:|------:|------:|------:|------------:|
| P-192                       |  192 | 0.250 | 0.768 | 0.911 | 1.000 | 0.647 |
| **P-224**                   |  224 | 0.200 | 1.000 | 0.875 | 1.000 | **0.647** |
| P-256                       |  256 | 0.500 | 0.850 | 0.848 | 1.000 | 0.775 |
| P-384                       |  384 | 0.300 | 0.947 | 0.966 | 1.000 | 0.724 |
| P-521                       |  521 | 0.000 | 0.992 | 0.964 | 1.000 | **0.001** |
| secp256k1                   |  256 | 0.667 | 0.922 | 0.906 | 1.000 | 0.864 |
| secp192k1                   |  192 | 1.000 | 0.969 | 0.870 | 1.000 | 0.958 |
| secp224k1                   |  224 | 1.000 | 0.814 | 0.987 | 1.000 | 0.947 |
| brainpoolP192r1             |  192 | 1.000 | 1.000 | 0.938 | 1.000 | 0.984 |
| brainpoolP224r1             |  224 | 1.000 | 0.911 | 0.955 | 1.000 | 0.966 |
| brainpoolP256r1             |  256 | 1.000 | 0.969 | 0.984 | 1.000 | **0.988** |
| brainpoolP384r1             |  384 | 1.000 | 0.911 | 0.940 | 1.000 | 0.962 |
| brainpoolP512r1             |  512 | 1.000 | 0.957 | 0.936 | 1.000 | 0.973 |
| FRP256v1                    |  256 | 1.000 | 0.856 | 0.961 | 1.000 | 0.952 |
| SM2                         |  256 | 0.500 | 1.000 | 0.941 | 1.000 | 0.828 |
| GOST-CryptoPro-A            |  256 | 0.667 | 0.883 | 0.918 | 1.000 | 0.857 |
| GOST-CryptoPro-B            |  256 | 0.833 | 1.000 | 0.996 | 1.000 | 0.955 |
| GOST-CryptoPro-C            |  256 | 1.000 | 0.921 | 0.934 | 1.000 | 0.963 |
| GOST-tc26-256-paramSetA     |  256 | 0.667 | 0.879 | 0.953 | 1.000 | 0.864 |
| GOST-tc26-512-paramSetA     |  512 | 0.214 | 1.000 | 0.963 | 1.000 | 0.674 |
| GOST-tc26-512-paramSetB     |  512 | 0.143 | 0.969 | 0.971 | 1.000 | 0.605 |

## Findings

- **The Nikolaev slide claim reproduces.**  P-224's special-prime
  score (0.200) is **the lowest among NIST curves** other than P-521
  (which is a true Mersenne prime, weight 2 — a different kind of
  structure).  Its overall PKM resistance (0.647) is **tied for
  lowest** with P-192.

- **Brainpool dominates the table.**  Every Brainpool curve scores
  ≥ 0.962 — confirming that the RFC 5639 pseudorandom-prime
  construction does what it was designed for.

- **The GOST tc26 512-bit curves are anomalous.**  Both
  paramSetA (0.674) and paramSetB (0.605) have very low Solinas weight
  for primes of their bit-length.  The Nikolaev criterion would flag
  them; whether this matters in practice depends on the (not-yet-
  understood) translation from "divisor sets exist" to "ECDLP
  feasible", which is the open question on the slide.

- **secp256k1 outperforms P-256.**  Bitcoin's curve scores 0.864 vs.
  P-256's 0.775 on the PKM axis.  Its prime `p = 2^256 − 2^32 − 977`
  has Solinas weight 6 vs. P-256's 5.  This does *not* mean secp256k1
  is more secure overall (P-256 has formal-derivation seeded
  parameters); it just means the specific PKM-attack handle is
  slightly worse on P-256.

## What this can't tell you

- **(1) is the only Nikolaev-specific signal.**  The other four are
  textbook ECDLP smoothness checks I bundled in for context.  A real
  PKM-criterion implementation would need the divisor-set search at
  cryptographic scale, which is currently infeasible.

- **The geometric mean is a UX choice.**  A curve with overall score
  0.99 could still be broken by a future PKM variant that exploits a
  signal we don't measure.  Use the column breakdown, not the
  aggregate, for decision-making.

- **P-521's 0.001 overall score is misleading.**  Mersenne primes
  enable *different* attack handles (e.g. fast modular reduction) than
  the PKM divisor-set attack targets — the two often live on different
  axes.  Don't read this as "P-521 is broken"; read it as "P-521's
  prime is structurally unique and deserves a separate analysis."

## Reproducing

```bash
cargo run --example pkm_audit_demo --release
cargo test --lib cryptanalysis::pkm_criterion
```

Tests include the **two formal claims**:

- `nikolaev_claim_p224_more_vulnerable_than_secp256k1` — asserts
  P-224's special-prime score is strictly below secp256k1's.
- `nist_curves_have_large_embedding_degree` — sanity-checks the MOV
  threshold check.

## References

- V. D. Nikolaev, *On the correctness of criterion for resistance to
  Petit-Kosters-Messeng method*, CTCrypt 2024.
- C. Petit, M. Kosters, A. Messeng, *Algebraic approaches for the
  elliptic curve discrete logarithm problem over prime fields*, PKC 2016.
- J. Solinas, *Generalized Mersenne numbers*, CACR-99-39, 1999.

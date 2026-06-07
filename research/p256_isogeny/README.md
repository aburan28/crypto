# P-256 Same-Order Isogenous Curve Search

Pure-Python (no Sage/PARI) implementation of the **first stages** of the
"Special Same-Order Isogenous Curves to P-256" search plan. Everything here
runs at full P-256 scale with the standard library only, and self-validates
its arithmetic against the published P-256 base point (`[n]G == O`) and against
brute-force point counting at toy scale.

This covers the recommended first implementation steps of the plan:
Phase 0 (invariants + classical rejection checks), Phase 1 (small-degree
isogeny-graph structure + shallow crawl), and the start of Phase 2 / Phase 4
(same-order model construction + cheap invariant scoring).

## Files

| file | purpose |
|------|---------|
| `ecfp.py` | field/curve arithmetic, polynomial root-finding mod p, Vélu 2- and odd-ℓ isogenies, division polynomials, order verification |
| `params.py` | canonical NIST P-256 domain parameters |
| `phase0_invariants.py` | Phase 0 — invariants + classical weak-curve rejection checks |
| `phase1_structure.py` | Phase 1 — exact small-prime isogeny-graph structure from Δπ, confirmed by direct construction for ℓ=2,3 |
| `phase1_crawl_ell2.py` | Phase 1 (Mode A/D) — BFS 2-isogeny crawl via Vélu, writes `curves.jsonl`/`edges.jsonl` |
| `phase1_crater_walk.py` | Phase 1 (Mode C/D) — odd-ℓ (Schoof-kernel) neighbour catalog + split-prime crater walks |
| `semaev.py` | Semaev summation polynomials over F_p (toy scale) |
| `phase7_relation_yield.py` | Phase 7.2 — Semaev relation-yield harness with controls (root vs same-order neighbour vs random) |
| `groebner.py` | minimal Buchberger/grevlex Gröbner engine over F_p (F4 degree-of-regularity proxy) |
| `phase7_groebner_dreg.py` | Phase 7.3 — solving-degree probe on Semaev decomposition systems, vs same-shape generic control |
| `phase5_endomorphism.py` | Phase 5 — endomorphism / GLV audit (CM norm-form reduction, automorphism + self-loop scan) |
| `analysis.py` | shared Phase 4 cheap-invariant scoring |
| `test_p256_isogeny.py` | self-tests, incl. brute-force validation of Vélu / Schoof codomains and the Semaev oracle |
| `*_output.txt` | captured reference outputs |

### Odd-ℓ isogeny construction (Schoof-style)

For odd ℓ the kernel of a rational ℓ-isogeny need not contain a rational point;
it is a Frobenius eigenspace in `E[ℓ]`. `ecfp.kernel_polynomials` finds, for each
eigenvalue λ (root of `X²−tX+p mod ℓ`), the kernel polynomial
`h_λ = gcd(xᵖ·D_λ − N_λ, ψ_ℓ)` of degree `(ℓ−1)/2`, and `velu_from_kernel_poly`
builds the codomain from the kernel's power sums `p₁,p₂,p₃` (all Vélu needs).
Validated against brute-force point counting for every ℓ∈{2,3,5,7,11,13} codomain
at toy scale.

## Run

```sh
python3 phase0_invariants.py        # invariants + verdict
python3 phase1_structure.py         # isogeny-graph structure + real 3-isogenous neighbour
python3 phase1_crawl_ell2.py 4      # 2-isogeny BFS to depth 4 -> curves.jsonl / edges.jsonl
python3 test_p256_isogeny.py        # all self-tests (<1s)
```

## Findings so far

### Phase 0 — classical weak-curve buckets are ruled out (as expected)

All four classical weaknesses are **same-order invariants**, so if P-256 lacks
them, no same-order neighbour over F_p can acquire them:

| check | result |
|-------|--------|
| anomalous (`n == p`) | OK (ruled out) |
| supersingular (`t ≡ 0 mod p`) | OK (ruled out) |
| smooth order (`n` composite) | OK — `n` is prime |
| small embedding degree (MOV/Frey-Rück, `k<100`) | OK — none found |

Recorded invariants: `t` (127-bit trace), `Δπ = t²−4p` (258-bit, negative),
small-prime valuations `v₃(Δπ)=v₅(Δπ)=1`, all others 0 over the search set.
Twist order is composite with small factors `3·5·13·179·(large prime)`.

### Phase 1 — the small-degree rational isogeny neighbourhood is *thin*

Because P-256 has **prime order**, it has no rational 2-torsion, so the rational
2-isogeny graph is a **single isolated vertex** (confirmed: the depth-4 Vélu
crawl finds 0 neighbours). More generally the rational ℓ-isogeny degree is fixed
exactly by `Δπ` (a same-order invariant) via the Deuring/Kohel volcano picture.
For `ℓ ∈ {2,…,31}`:

| class | primes ℓ | rational ℓ-isogenies |
|-------|----------|----------------------|
| **ramified** (`ℓ \| Δπ`) | 3, 5 | 1 each (vertical) |
| **split** (`(Δπ/ℓ)=+1`) | 11, 13, 17, 23, 29 | 2 each (horizontal crater) |
| **inert** (`(Δπ/ℓ)=−1`) | 2, 7, 19, 31 | 0 |

All volcano heights are 0 over this set (`v_ℓ(Δπ) ≤ 1`), i.e. P-256 sits on the
crater for every small ℓ — there is essentially no vertical depth to exploit.
The ℓ=2 and ℓ=3 predictions were confirmed by direct construction (2-division
roots; ψ₃ roots + Vélu).

### First constructed same-order neighbour (3-isogenous)

A real, order-verified same-order neighbour of P-256:

```
a = 69813713139673779928897632679399297923312904910852909497923264316660085700227
b = 46459849483693785433638644630303853891775023127090599264446658487217616473547
j = 60359795834994757875819835712620686501669024665003573244969414519504858434740
#E'(F_p) = n  (verified by random-point order test)
```

Its **cheap invariant score is zero**: `j ∉ {0,1728}`, `a ≠ −3`, dense
`a,b` (Hamming weight 121/136), automorphism group of size 2 — a completely
generic ordinary model. This is the H0 ("generic hardness") baseline: the first
non-trivial neighbour is indistinguishable from a random ordinary curve, exactly
as the plan anticipates for the overwhelming majority of nodes.

### Split-prime crater walk (Phase 1 Mode C/D)

`phase1_crater_walk.py` constructs real same-order neighbours via odd-ℓ isogenies
and walks the horizontal crater. Results:

* **12-neighbour catalog** (one step per ramified/split ℓ ≤ 29): every codomain
  is order-verified `= n`, distinct from the root, and **cheap-score zero**
  (no special `j`, no `a=−3`, dense coefficients, automorphism group 2). The
  per-ℓ neighbour count exactly matches the Δπ predictor — an independent
  Schoof-eigenvalue confirmation (5→1, 7→0, 11/13/17/23/29→2).
* **Horizontal walks** of 7 (ℓ=11) and 5 (ℓ=13) distinct same-order nodes, all
  generic.
* **Finding:** P-256's defining sparse model `a=−3` is a property of the *root
  only* — none of its isogenous neighbours inherit it. The "special" sparse model
  is not isogeny-invariant.

### Phase 7 — Semaev relation yield: no signal (H0 confirmed)

`phase7_relation_yield.py` reduces P-256 mod a toy prime `q`, builds its rational
same-order ℓ-isogenous neighbour, and compares **index-calculus relation yield**
(fraction of random targets decomposable into `m` factor-base points) against
strict controls, with an *identical* factor-base size across all curves
(same-work comparison) and a Semaev-S₃ oracle cross-check:

| q | root | neighbour | same-order ctrl (mean) | random ctrl (mean) | neighbour/root |
|---|------|-----------|------------------------|--------------------|----------------|
| 12007 | 0.245 | 0.263 | 0.266 | 0.276 | 1.07 |
| 40009 | 0.105 | 0.083 | 0.089 | 0.088 | 0.79 |

The neighbour's relation yield is **statistically indistinguishable** from the
root and from same-order / random controls (ratios bracket 1.0, within control
spread), and the signal does not grow with field size. This is the **H0
"generic hardness"** outcome: isogenous movement does **not** change the
relation-generation geometry — no Level-1 signal (would require ≥5×).

### Phase 7.3 — Gröbner degree of regularity: no anomaly (H3 absent)

`phase7_groebner_dreg.py` runs a pure-Python Buchberger (grevlex) on the
2-point factor-base decomposition system
`{ S₃(x₁,x₂,x_R)=0, g(x₁)=0, g(x₂)=0 }` (g = factor-base vanishing polynomial)
and records the **solving degree** (max remainder/S-poly degree), reduction count
and basis size. Over q=12007, |FB|=5, 10 consistent targets per curve:

| curve | Semaev solving degree | reductions | generic same-shape control |
|-------|-----------------------|------------|----------------------------|
| root (P-256 mod q) | 9.0 | ~2095 | 9.0 |
| same-order neighbour | 9.0 | ~2122 | 9.0 |
| 3× same-order control | 9.0 | ~2104 | 9.0 |
| 3× random control | 9.0 | ~2131 | 9.0 |

Every curve — including a **generic polynomial of identical monomial support** —
solves at the *same* degree with the *same* work. The solving degree is fixed by
the system shape, not the curve: no curve-dependent degree-of-regularity drop, no
extra low-degree syzygy. H3 (Semaev/Gröbner anomaly) is absent.

### Phase 5 — endomorphism / GLV audit: no usable endomorphism (H1 absent)

`phase5_endomorphism.py` reduces the CM norm form `[1, t, p]` (norm of
`a+b·π`). The reduced form is `[1, 1, ≈2²⁵⁶]`: the smallest **non-scalar**
endomorphism has degree ≈ 2²⁵⁶. The only efficiently computable endomorphism is
Frobenius (degree p), which acts as the **identity** on the rational order-n
subgroup (eigenvalue λ=1) — useless for GLV. All 12 catalogued neighbours have
automorphism group size 2, none has special j∈{0,1728}, and none is a self-loop.

→ **rho bits saved ≈ 0** (below the "interesting" threshold of 1 bit). Unlike
secp256k1 (j=0, CM disc −3, genuine GLV), P-256 and its same-order neighbours
have no GLV-style endomorphism.

## Correctness

`test_p256_isogeny.py` (8 tests, <1 s) independently validates:
* Vélu 2-/3-isogeny codomains vs brute-force point counting (236 toy codomains);
* odd-ℓ Schoof-kernel isogenies ℓ∈{5,7,11,13} vs brute force (106 codomains),
  incl. kernel-polynomial degree `(ℓ−1)/2`;
* the Semaev S₃ oracle (vanishes on real collinear triples);
* the Buchberger engine (ideal membership) and binary-quadratic-form reduction;
* the P-256 Phase 0 / Phase 1 facts (`[n]G=O`, 0 rational 2-isogenies, 1 verified
  same-order 3-isogeny).

## Status vs the plan

Done: Phase 0; Phase 1 (structure + ℓ=2 crawl + odd-ℓ catalog + crater walks);
Phase 2 (same-order model construction + order verification); Phase 4 (cheap
scoring); **Phase 5 (endomorphism/GLV audit)**; Phase 7.2 (relation yield) and
**Phase 7.3 (Gröbner solving degree)** with Phase 3 controls and Phase 13 red-team
checks (same-work, same-shape generic control, controls, scaling).

Not yet implemented:
* Phase 6 deep conductor audit (full factorisation of Δπ beyond small primes);
* Phase 7.4 Betti/syzygy tables; Phase 8 low-genus correspondence search;
* Phase 9 explicit isogeny-map transport; Phase 10 scoring sigmoid roll-up.

**Bottom line so far:** every measurable probe — classical invariants, volcano
structure, cheap invariants, toy-scale relation yield, Gröbner solving degree, and
the CM endomorphism lattice — places P-256's same-order isogenous neighbours
squarely in the generic-ordinary bucket. The vertical volcano is trivial (height 0
at all small ℓ), the horizontal crater nodes are generic, no GLV endomorphism
exists, and isogenous movement leaves both relation yield and degree of regularity
unchanged. No hypothesis H1–H4 shows any signal; H0 (generic hardness) holds
across the board.

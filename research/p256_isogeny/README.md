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
| `test_p256_isogeny.py` | self-tests, incl. brute-force validation of Vélu codomains |
| `*_output.txt` | captured reference outputs |

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

## Correctness

`test_p256_isogeny.py` validates the Vélu codomain formulas independently by
brute-force point counting over 236 toy codomains (every 2- and 3-isogeny found
for small `a,b` over `p ∈ {101,1009,2003}`), confirming order preservation, and
re-checks the P-256 Phase 0/Phase 1 facts. All pass in <1 s.

## Not yet implemented (next steps from the plan)

* odd-ℓ isogeny construction for ℓ ≥ 5 (kernel-polynomial factoring of ψ_ℓ, or
  modular-polynomial roots) to walk the split-prime crater cycles (11,13,17,23,29);
* Phase 3 random/same-order controls;
* Phase 5 endomorphism/GLV audit; Phase 6 conductor audit;
* Phase 7 Semaev/Gröbner harness at toy + reduced-modulus scale (highest value);
* Phase 8 low-genus correspondence search.

The structural conclusion driving prioritisation: the **vertical** volcano
structure around P-256 is trivial at small primes (height 0 everywhere), so any
non-generic behaviour would have to come from the **horizontal** crater nodes or
from the algebraic (Semaev/Gröbner) geometry — not from conductor depth.

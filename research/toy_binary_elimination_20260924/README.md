# Tiny isogeny-class and Boolean elimination audit

Dependency: [PR #708](https://github.com/aburan28/crypto/pull/708).
Classification: accounting; stage diagnostic. No speedup or regularity claim.

## Results

We enumerated 53 curves y^2+xy=x^3+c over GF(8), GF(16), and GF(32), using
all nonzero c and a2=0. This covers 13 point-count classes within this family,
not all twists or all binary elliptic curves. Every affine point is recorded.
Equal point counts certify membership in the same finite-field isogeny class;
no explicit isogeny maps, degrees, or volcano strata are constructed.

For each curve we tested V=span(1,t) and V squared, with six Boolean variables
encoding ordered triples. All 53 Frobenius transports preserved exact zero
and rational-lift assignment sets. All 40,064 point-pair checks passed.

| Field | Systems | Unsatisfiable systems | Linear-consequence completion degree |
|---|---:|---:|---|
| GF(8) | 14 | 0 | 5 |
| GF(16) | 30 | 10 | 4–5 |
| GF(32) | 62 | 36 | 4–5 |

Across the 106 systems there were 306 polynomial-zero triples and 206 triples
admitting rational zero-sum lifts. These totals include the transported control
systems and repeated/zero coordinates; they are not independent relation yields.
Of 46 unsatisfiable systems, contradiction certificates appeared at D=3 in 8,
D=4 in 26, and D=5 in 12. The remaining 60 systems had no contradiction.
Linear consequences were complete at D=4 in 40 systems and D=5 in 66.

The degree here is the explicitly specified Boolean Macaulay filtration
threshold. It is neither a measured F4/F5 solving degree nor a proved degree
of regularity. Fixed row ordering and generator encoding affect the diagnostic.

## A misleading comparison caught by rational lifting

In GF(32) with modulus t^5+t^2+1, c=24 and c=30 give non-Frobenius-conjugate
curves with 36 rational points. Both fixed-subspace systems have one polynomial
zero and linear-completion degree 4. Yet only c=24 admits a rational lift.

| c | Ordered x triple | Polynomial zeros | Rational-lifted triples | Matrix row XORs through D=9 |
|---|---|---:|---:|---:|
| 24 | (2,2,2) | 1 | 1 | 1815 |
| 30 | (3,3,3) | 1 | 0 | 2725 |

These are the only fixed-basis, nonconjugate, equal-point-count,
equal-positive-polynomial-zero-count pair in this panel. The remaining nine
equal-zero-count pairs are unsatisfiable. Thus this panel contains **no such
pair with matching positive rational-lift counts**. A row-XOR difference here
cannot support a claim of cheaper verified relation generation. The repeated
point solution is also unsuitable as evidence of generic decomposition yield.

## Why the original saturation metric was refined

The first run reached full ideal saturation at D=6 in every case. This has an
unavoidable ambient-dimension floor: every nonzero ideal of the Boolean function
algebra contains a point indicator, whose unique squarefree polynomial has
degree six. Rows admitted at D<6 cannot span that indicator. An observed common
threshold of six therefore does not measure equal Groebner difficulty.

The first protocol, source, and results are retained. FOLLOWUP_PROTOCOL.md was
written before the refined run, which additionally records contradiction and
linear-consequence thresholds. All fields and counters from the original run
match the refined run exactly. Projection-rank diagnostics are excluded from
the original row-XOR count, which measures only elimination, not total cost.

## Verification and replay

From the repository root:

```sh
python3 research/toy_binary_elimination_20260924/verify.py
```

The verifier replays both stages, checks all original fields are preserved,
independently derives the dimension of linear consequences from the exhaustive
solution sets using dense binary linear algebra, checks three elementary ideals,
and validates source/input/result hashes. Full Boolean ideal rank is independently
bounded by exhaustive zero counts; generated rows are checked on every zero.

The first verifier attempt compared Python point tuples directly to JSON lists
and failed at replay equality. Normalizing the computed object through JSON
fixed that representation-only verification bug. Neither frozen result changed.

## Decision

Retain this as a correctness and metric-selection control. No measured
isogeny-descent advantage has been established. Explicit neighboring isogenies,
nondegenerate rational decompositions, matched verified workloads, and actual
F4/F5 traces remain outstanding. This bounded audit should not be extrapolated
to ECC2K-130. All end-to-end costs, rho ratios, speedups, and regularity fields
remain null.

Background on finite-field isogeny classification: [Jordan et al., Abelian
varieties isogenous to a power of an elliptic curve](https://math.mit.edu/~poonen/papers/power-of-E.pdf),
citing Tate's isogeny theorem. For an elliptic curve the characteristic
polynomial is X^2-(q+1-#E)X+q, so the point count fixes it.

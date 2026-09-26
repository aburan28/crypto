# Frobenius quotient follow-up: toy decompositions and a benchmark correction

Date: 2026-09-26 UTC. Stage diagnostic only: end-to-end cost, `S`, and speedup
are **unset**. No work below `2^61` has been demonstrated.

This directory records a follow-up note on target-only four-point decompositions
on `E: y² + xy = x³ + 1` over `F_{2^n}`. It also independently checks every
claim in that note that can be checked without the note's solver code.

```sh
python3 research/frobenius_quotient_followup_20260926/verify_claims.py > /tmp/out.json
```

Pure Python 3, no dependencies, seeded (`20260926`). It takes about 3 minutes on
one core. The committed output is `results.json`.

## Provenance and status

The SAT measurements in the table under **Reported solver runs** were supplied
as a written note. **The CryptoMiniSat/Z3 models, raw solver logs and target
files are not in this repository.** Those numbers are recorded as reported.
They are not reproduced here, and no downstream claim may treat them as
reviewed evidence until the source and raw runs are committed. The algebraic
identities and counting claims *are* re-derived here by independent code; see
**Independent verification**.

## Boundaries (stated before interpreting)

- **Floor (counting, not heuristic).** Fix four distinct point orbits under
  sign and Frobenius, each of size `2n`, in a subgroup `H`. For a target `T`
  uniform on `H`, let `N(T)` count the assignments (one point per orbit) that
  sum to `T`. Then `E_T[N(T)] = (2n)^4/|H|` and `Pr_T[N(T) > 0] ≤ (2n)^4/|H|`.
  These are exact. They cover fixed labels and a uniform target, not labels
  chosen using the target.
- **Reference.** Pollard rho on the ECC2K-130 odd subgroup (`|H| ≈ 2^129`).
  None of the runs below approaches a comparison with it.

## Formulations (as reported)

- **Line formulation.** Write `R = (a, b)`, and let `t`, `u` be the
  x-coordinates of the intermediate sums `A`, `B`, with `A`, `B`, `−R`
  collinear (so `A + B = R`). Let the line have slope `ℓ` and intercept
  `ν = b + a + aℓ`. Then

  `u = t + ℓ² + ℓ + a`, and `t·u = aℓ² + a² + b + a`.

  After substitution this is one quadratic Boolean relation in `(t, ℓ)`. It is
  a geometric reparameterization that avoids the degree-14 quotient inverse. It
  does not lower the Boolean degree below the target-specialized `f3(t,u,a)`
  (already quadratic once `a` is fixed), and it does not remove phases or
  representatives. No novelty is claimed.
- **Quadratic subgroup certificate.** For `x ≠ 0` with `Tr(x) = 0`, add an
  auxiliary `q` with `q² + √x·q + 1 = 0` and `Tr(q) = 0`. This gives
  `x = (q + 1/q)² = x(2Q)` with `Q ∈ 2E`, so `x` is the x-coordinate of a point
  in `4E`. Conversely, every point of `4E` has such a `q`. This moves the
  subgroup rejection into the search. It does not encode a key list.

## Reported solver runs (not reproduced here)

These are CryptoMiniSat 5.16.0 runs on one thread with a 2 s solve budget
(model construction excluded from the budget but included in the totals).
All four x-coordinates are distinct. Only the target is supplied; every
returned tuple is checked on the curve. Each cell reads *verified solutions /
targets; total seconds*.

| Target set     | Formulation     | n=13, s=4     | n=19, s=6     |
|----------------|-----------------|---------------|---------------|
| Planted        | specialized_xor | 0/4; 8.083 s  | 0/4; 8.158 s  |
| Planted        | specialized_cnf | 0/4; 8.170 s  | 0/4; 8.066 s  |
| Planted        | subgroup_xor    | 1/4; 6.997 s  | 0/4; 8.165 s  |
| Planted        | line_xor        | 1/4; 6.889 s  | 0/4; 8.102 s  |
| Planted        | line_subgroup   | 1/4; 7.879 s  | 0/4; 8.133 s  |
| Earlier random | line_xor        | 0/4; 8.051 s  | 0/4; 8.131 s  |
| Earlier random | line_subgroup   | 1/4; 7.779 s  | 0/4; 8.142 s  |

That is 4 of 56 unknown-phase runs verified (three distinct degree-13
targets) and 0 of 28 at degree 19. There is one timing trial per cell, so these
runs establish no speedup. Native XOR against expanded CNF showed no difference
at this budget. In 16 further fixed-phase diagnostics, 1 solved; they are
excluded from the target-only counts.

## Independent verification (this directory)

All figures below come from `results.json`.

| Claim in the note | Independent result |
|---|---|
| Odd subgroup order for n=131 is 680564733841876926932320129493409985129 | matches (`#E = 4·h`, `h` probable prime, `log2 h ≈ 129.0`) |
| n=29 odd subgroup `4E` has order 134207651, composite | matches (`#E = 536830604`) |
| `E_T[N] ≈ 2^−96.8663` at n=131 with `262^4` assignments | `log2 = −96.86631` |
| Certificate accepts exactly 1001 of 8191 nonzero x at n=13, matching `4E` | 1001 accepted, 8191/8191 agree with true membership (`|4E| = 2003`) |
| Certificate matches membership (n=19) | 4000/4000 random x agree (519 accepted) |
| Line identities hold, including tangents, n=13 and 19 | 100/100 each (20 tangent cases each); the points are reconstructed from `(t, ℓ, ν)` |
| Exact means 228.14578 / 15.93251 / 0.08432 | 228.14578 / 15.93251 / 0.08432 (`|H|` = 2003 / 130873 / 134207651) |

The density check uses an independently drawn random factor space per field,
with the basis recorded in `results.json`, and exhaustive pair matching over 8
targets per tuple:

| n, s | valid orbits | tuples | positive cases | observed mean | uniform-target mean |
|------|---:|---:|---:|---:|---:|
| 13, 4 | 4  | 1  | 8/8     | 231.375 | 228.146 |
| 19, 6 | 10 | 32 | 256/256 | 15.520  | 15.933  |
| 29, 7 | 16 | 16 | 10/128  | 0.0781  | 0.0843  |

This reproduces the note's qualitative correction on a different factor space:
compatible phase assignments saturate at n=13 and n=19 and become sparse at
n=29. Grid cells are correlated, not independent trials.

## Consequence for earlier benchmarks

For n=13 and n=19, any four distinct valid orbit labels almost always admit a
lift to a given target. So success at those sizes does not test the
target-dependent selection of compatible keys, which is the hard part at full
size (`≈ 2^−96.9` expected lifts per arbitrary label tuple at n=131). n=29 is
the smallest size in this set that gives an informative sparse control. No
degree-29 SAT solve is claimed. The counting bound gives no lower bound against
an adaptive algebraic solver.

## Not measured

Relation-matrix rank, Groebner solving degree, repeated timing trials, and
end-to-end discrete-log cost were not measured.

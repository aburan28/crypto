# Radical pair-image results — engineering diagnostic, 2026-09-21

All 276 accepted msolve processes complete, and all 34 target certificates pass.
The radical presentation lowers maximum F4 degree on every target.  This is a
solver-stage result; no complete-ECDLP cost or attack speedup is established.

## Primary end-to-end metric

| Variant | Class | Total common operations | S | Matched rho ratio | Generic floor ratio | Correctness |
| --- | --- | --- | --- | --- | --- | --- |
| Direct sparse S3 tree | engineering reference | null | null | null | null | all certificates pass |
| Two norm equations | engineering control | null | null | null | null | all prime certificates pass |
| Radical pair image | engineering candidate | null | null | null | null | all certificates pass |

Null means unmeasured.  Relation collection, final relation-matrix work,
individual-log extraction and a common-operation conversion are absent.

## Algebraic boundary

Let `W` contain `b` distinct field elements and

```text
P_W = {(x+y,xy): x,y in W}.
```

There are `b(b+1)/2` points and the same number of monomials `s^i t^j`
with `i+j <= b-1`.  If a polynomial of degree at most `b-1` vanishes on
`P_W`, substituting `s=x+y,t=xy` gives a polynomial of degree below `b`
separately in `x,y` that vanishes on `W^2`.  Bivariate interpolation makes it
zero.  Therefore those monomials form the quotient basis and the radical image
ideal has a DRL basis of `b+1` degree-`b` relations:

```text
H(T) = 1 + 2T + ... + b T^(b-1),       d_reg(radical image) = b.
```

Ordered membership has top ideal `<x^b,y^b>` and regularity `2b-1`.

| Presentation | Membership regularity | Ratio to exact image floor b |
| --- | ---: | ---: |
| Ordered `h(x)=h(y)=0` | `2b-1` | `(2b-1)/b` |
| Radical unordered-pair image | `b` | `1.000` |

The radical basis is generated once by eliminating `x,y` from
`<h(x),h(y),s-x-y,t-xy>` and converting to DRL.  The two norm equations have
the right set of roots but retain ordered/diagonal scheme multiplicity; they do
not consistently improve the solver.

## Matched F4 stage

The stage table includes producing targets only.  `Wall C/R` is the geometric
mean of paired radical/direct process-time ratios over three repetitions per
target; the interval is a 4,000-draw paired bootstrap.  Matrix ratio is
candidate/reference for the median peak F4 round area.  `Cold batch C/R`
charges radical preprocessing once and includes all targets in the cell.

| Cell | b | Producing / targets | F4 degree R→C | Matrix C/R | Wall C/R [95% CI] | Cold batch C/R |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| prime-main-b3 | 3 | 2/3 | 6→5 | 1.0111 | 0.7095 [0.4253, 0.9777] | **1.9189** |
| prime-main-b5 | 5 | 2/3 | 9→7 | 0.1398 | 0.3737 [0.3525, 0.3962] | 0.4958 |
| prime-main-b7 | 7 | 2/3 | 11→9 | 0.2318 | 0.2062 [0.1982, 0.2133] | 0.2386 |
| prime-main-b9 | 9 | 2/3 | 15→11 | 0.1942 | 0.1449 [0.1254, 0.1633] | 0.1490 |
| prime-holdout-b5 | 5 | 2/3 | 9→7 | 0.0551 | 0.4750 [0.4237, 0.5434] | 0.6414 |
| prime-holdout-b7 | 7 | 2/3 | 11→9 | 0.5807 | 0.2043 [0.2016, 0.2071] | 0.2507 |
| prime-holdout-b9 | 9 | 2/3 | 13→11 | 0.3077 | 0.1173 [0.1115, 0.1238] | 0.1327 |
| prime-holdout-b11 | 11 | 2/3 | 15→13 | 0.3116 | 0.1083 [0.0841, 0.1418] | 0.1399 |
| binary-k2 | 3 | 3/3 | 5→4 | 0.2508 | 0.8322 [0.7523, 0.9304] | 0.9488 |
| binary-k3 | 7 | 7/7 | 6,7→5 | 0.0345 | 0.3099 [0.2658, 0.3561] | 0.3130 |

The `b=3` prime cell is an explicit regression once setup is charged.  It is not
promoted.  Every `b>=5` prime main/holdout cell and both binary cells have cold
batch ratio below one in this solver-stage suite.  Timing remains secondary:
shared-host variance is visible in raw repetitions, especially at `b=11`.

## Correctness and scope

- Direct and symmetric finite root sets agree on every target after quotient
  and exact quadratic-root recovery.
- Recovered roots belong to the identical factor base.
- Finite chain roots equal every nondegenerate elliptic-curve group relation.
  Omitted `P,-P` intermediate-infinity branches reduce to lower-arity relations.
- Prime main uses `GF(17), y^2=x^3+2x+2`; holdout uses
  `GF(29), y^2=x^3+x+1`.
- Binary controls use `GF(2^4)/GF(2^2)` and `GF(2^6)/GF(2^3)` with the
  nonzero subfield factor base, excluding the ordinary curve's `x=0` two-torsion.
- The prior `GF(2^8)/GF(2^4)` feasibility run is censored after 90 seconds.
  No direct comparison was launched and no degree/runtime claim is inferred.

Accepted run 002 stores 276 raw process logs, 34 certificates, every msolve
input, source/input/log hashes and full F4 rounds.  Independent audit status is
`pass`. Run 001 is superseded because it stored absolute raw-log paths. Run 000
remains rejected because `msolve -g 0` measured full solution parametrization
rather than the intended Gröbner stage.

## Interpretation

This round establishes an exact presentation improvement and finite F4-stage
engineering gain.  It does not establish a relation-yield improvement, a
complete index-calculus gain, an exponent change, a rho crossover or deployed-
curve security impact.  The next required experiment is a cold complete-DLP
pipeline using this presentation, with relation collection, rank, final linear
algebra and individual-log extraction priced in a common operation unit.

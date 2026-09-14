# Pulling the coefficient equation back to support image space

The new Nagao solver rejects impossible function coefficients **before**
generating them. It preserves every admissible coefficient set on the exhaustive
five-bit panel. The matched campaign contains 96 solver cells and 18 complete
tiny ECDLP runs. All completed answers passed independent checks.

This is an **engineering candidate, classification pending normalized cost**.
There is no exponent improvement or established end-to-end speedup. The tiny
rho reference remains cheaper in every recorded field-operation component.
This is an independently derived implementation idea, not a literature-novelty
claim. It applies to the binary a=0 chart below; prime fields, a=1 Koblitz curves,
and F4/Gröbner engines have not been changed or benchmarked in this round.

## Frozen boundary and decision

[contract.json](contract.json) predates the stage campaign. The baseline is
`solver_08/image_solver.py`; the controls are the existing symmetric S4 and
chained S3 CryptoMiniSat adapters. The WDSat regression uses a=1 and another
polynomial protocol, so this freezes an equivalent matched a=0 suite under the
[parent accounting contract](../../index_calculus_baseline_20260914/ec_index_calculus_contract.json).
It does not replace the accepted WDSat baseline.

There remain O(2^(2d)) conditioned (h2,z) branches, with at most two b values
per h2 and two image preimages per branch. The dense binary elimination costs
O(d²) counted field additions and O(d) field multiplications per branch, giving
an O(d² 2^(2d)) field-API upper bound for this implementation, apart from
field-degree-dependent arithmetic and extraction. No exponent is fitted from
three sizes or censored cells.

For B signed factor points, a uniform affine target can be represented by at
most 8*C(B/2,3) signed triples of distinct abscissas. Thus its success probability
is at most min(1, 8*C(B/2,3)/(#E-1)); domain exclusions only reduce it. This
counting ceiling is unchanged. It does not price the group operations or
linear algebra and is not a calibrated full-DLP cost floor. Planted supported
targets cannot estimate natural relation yield.

The frozen diagnostic target was 20% fewer counted field-API calls on complete,
equal-output pairs. Every one of the eight eligible stage pairs met that sum
target. **Multiplications increased** on the eleven-bit stage panel, however;
an unweighted API sum does not establish a calibrated total-cost gain. The
larger panels do not provide complete paired enumeration. Keep this as an
experimental solver; do not promote it as an attack improvement.

| Variant/reference | Full-DLP S | Cost/rho | Cost/floor | Correctness evidence | Classification |
|---|---|---|---|---|---|
| Quadratic image baseline | null | null | null | 6/6 tiny scalar recoveries | Reference |
| Coefficient pullback | null | null | null | 6/6 tiny scalar recoveries; exhaustive support check | Engineering candidate; pending calibration |
| Symmetric S4 | null | null | null | 43 exhaustive toy targets; 24 campaign cells | Stage control |
| Chained S3 | null | null | null | 43 exhaustive toy targets; 24 campaign cells | Stage control |
| Matched rho | null | null | null | 6/6 tiny scalar recoveries | Reference; common unit uncalibrated |

## Algebra implemented

Use E: y²+xy=x³+1, target (r,s), and the x-coordinate subspace V of dimension d.
The function is x²+a*x+c+b*y with c=r²+a*r+b*(r+s), so its norm contains
the target factor (x+r). Its monic residual cubic H has

```
h2 = b²+b+r
h1 = a²+a*b+r*h2
```

As in the baseline, choose h2 in V, solve for nonzero b, and condition on
z in V with z != 0,r,h2. Put t=r+z, u=h2+z, K=H_(a=0)(z), and
c0=r*h2+z*u. The two remaining roots lie in V precisely when their product
v lies in I_u = image(w -> w²+u*w, w in V). Because u is a nonzero member of V,
the kernel is {0,u} and dim(I_u)=d-1. The two coefficient conditions are

```
H_a(z) = K+t*a²+b*z*a = 0
v = a²+b*a+c0
```

Adding t times the second equation eliminates a²:

```
b*r*a = t*v+t*c0+K
gamma = t/(b*r)
delta = (t*c0+K)/(b*r)
a = gamma*v+delta
```

Substitution yields the binary-linear equation

```
C2*v²+C1*v = rhs, restricted to v in I_u
C2 = t*gamma²
C1 = b*z*gamma
rhs = K+t*delta²+b*z*delta
```

`imagePreimages` solves this restriction directly on a basis of I_u, then
produces a. C2 is nonzero, so the kernel dimension is at most one. All field
column and preimage additions are charged. Precomputed image squares and
inverses are charged to setup. The r=0 chart falls back to the original exact
generator. Inherited recovery rejects zero/repeated/target abscissas, checks
actual curve points, and verifies the signed sum. This still guesses one root;
it is a hybrid coefficient search, not a root-free solver.

## Measurements

Python 3.12.14, pycryptosat 5.14.7, one process at a time. Each stage cell has a
cooperative three-second budget including setup; overruns remain in raw data.
There is one repetition: times are descriptive, with no confidence-interval
or runtime-speedup claim. Two uniform and two known-decomposable fresh targets
per size exclude historical solver_02 through solver_09 targets. Those targets
are now development data and should not be reused as fresh holdouts.

The complete eleven-bit enumeration panel has four targets and seven unique
relations per variant. Counts include setup, search, extraction and verification.

| Variant | Field additions | Field multiplications | Field squarings | Field-API sum |
|---|---:|---:|---:|---:|
| Quadratic image | 142,960 | 86,086 | 72,461 | 301,507 |
| Coefficient pullback | 97,804 | 90,822 | 11,545 | 200,171 |

The API sum fell 33.61%, but multiplications rose 5.50%. On the exhaustive
five-bit validation, generated candidate occurrences fell from 2,798 to 444
with identical valid (a,b) sets and complete projected solution sets.

| (field bits, support dimension) | Variant | Resolved cells / 8 | Timed out / 8 |
|---|---|---:|---:|
| (11,5) | Each of all four variants | 8 | 0 |
| (23,7) | Coefficient pullback | 2 | 6 |
| (23,7) | Each other variant | 0 | 8 |
| (29,8) | Coefficient pullback | 1 | 7 |
| (29,8) | Each other variant | 0 | 8 |

The three additional completions are first relations on supported targets at
2.5664, 1.8957 and 2.8386 seconds. All large uniform cells and all large
enumerations timed out. Partial enumerations may contain valid relations;
they remain incomplete and are excluded from equal-output cost ratios.
These data establish neither a uniform-target advantage at those sizes nor
an asymptotic crossover. S3/S4 are SAT controls, not F4 measurements.

## Full cold ECDLP smoke tests

[e2e_contract.json](e2e_contract.json) fixes two tiny fields, three seeds, and
the two Nagao variants plus rho (18 runs). Every run rebuilds the curve/order
check, selects a full-order A, constructs G=4A and a planted Q=kG, builds the
factor base, collects independent relations, solves the modular matrix,
performs individual descent and verifies [k]G=Q. The secret is used only for
target setup and final verification. Exhaustive curve-order setup is charged
and intentionally makes this a correctness smoke test, not a scaling model.

Raw factor points need not belong to the prime subgroup. The harness projects
the **whole relation** by cofactor 4. For R=sA, the projected matrix RHS is
s mod N relative to G. Equal or opposite projected columns are deduplicated;
projected identities have known zero log. For individual R=Q+sA, the recovered
projected log L satisfies L=4k+s, so k=(L-s)/4 mod N. This avoids the incorrect
assumption that the normal-basis support is a prime-subgroup factor base.

All six paired IC paths used identical generated targets and first relations.
There are 4 and 10 matrix columns respectively; all reach full rank. A separate
post-run audit checks all 126 IC oracle calls, including empty answers, against
the independent group pair oracle (61 distinct targets). It also verifies
all recovered scalars independently. Audit work is validation overhead, not
an oracle available to the timed solver.

The following sums cover three cold runs per row. Field API components have
no measured common-operation conversion; scalar modular operations are
reported separately below. Neither sum is normalized S.

| Field bits / subgroup N | Variant | Field additions | Field multiplications | Field squarings | Cold seconds (secondary) |
|---|---|---:|---:|---:|---:|
| 5 / 11 | Quadratic image | 27,206 | 27,606 | 9,539 | 0.068100 |
| 5 / 11 | Pullback | 16,066 | 22,476 | 4,595 | 0.052301 |
| 5 / 11 | Rho | 2,037 | 3,996 | 717 | 0.008804 |
| 9 / 127 | Quadratic image | 226,874 | 199,396 | 105,191 | 0.765548 |
| 9 / 127 | Pullback | 115,830 | 153,397 | 39,495 | 0.507000 |
| 9 / 127 | Rho | 25,866 | 45,389 | 19,869 | 0.159190 |

| Field bits | Variant | Scalar modular additions | Multiplications | Inversion calls |
|---|---|---:|---:|---:|
| 5 | Each IC variant | 162 | 168 | 15 |
| 5 | Rho | 78 | 3 | 3 |
| 9 | Each IC variant | 762 | 993 | 33 |
| 9 | Rho | 141 | 27 | 3 |

Phase times, failed attempts, per-attempt counters, exact targets and scalar
certificates are in [e2e_results/raw.jsonl](e2e_results/raw.jsonl).
Field inversions are expanded by the inherited counter; scalar inversions are
native calls, not expanded. Python control flow, bit tests and coordinate
conversions are not represented by the field API vector. Cold wall time covers
those costs, but cannot substitute for the repository's calibrated unit.
There is no warm amortization claim, calibrated S, floor ratio or end-to-end
speedup. Rho uses a three-partition Floyd walk on the same subgroup and target.

## Reproduce and next gate

From the repository root, create a Python 3.12 environment with
`pip install pycryptosat==5.14.7`. Use new output names on every iteration;
the runners refuse to overwrite directories and the comparison refuses to
overwrite its output. No Magma or proprietary Gröbner solver is required.

```sh
python research/nagao_relations/coefficient_pullback/test_pullback.py -v
python research/nagao_relations/coefficient_pullback/experiment.py --output /tmp/pullback-stage-new
python research/nagao_relations/coefficient_pullback/e2e.py --output /tmp/pullback-e2e-new
python research/nagao_relations/coefficient_pullback/compare.py --stage /tmp/pullback-stage-new --e2e /tmp/pullback-e2e-new --output /tmp/pullback-comparison-new.json
```

The stage runner repeats the exhaustive 43-target check and 86 Semaev-control
checks before timing. Unit tests independently enumerate fibers of the
restricted linear map and check modular rank, dependent rows and inconsistent
relations. [comparison.json](comparison.json) freezes the post-run audit and
pair diagnostics. Source hashes in both run manifests identify measured code;
the stage manifest also pins the inherited solver dependencies. The later
comparison and unit tests are validation code, not measured solver changes.

Next: reduce the per-branch image elimination cost, preserving this exact
support test and counting its setup. Before promotion, freeze fresh uniform
holdouts with longer budgets, rerun the paired cold pipeline at larger groups,
calibrate field/scalar/control costs, and require the AGENTS.md normalized
total-cost gate. F4 comparison needs a separately pinned compatible backend;
the current SAT adapter cannot be relabelled as a Gröbner benchmark.

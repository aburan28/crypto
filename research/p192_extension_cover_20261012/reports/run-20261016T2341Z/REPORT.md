# P-192 prime-field cubic-cut optimization screen

**Result.** Two alternative cubic-cut implementations reproduced the exact
fixed-origin input data for 16 frozen two-pivot quotient classes, comprising
32 positive and negative residual cuts. Expanding extension-field cubic
equations into prime-field coefficient rows and computing the resultant
over \(\mathbb F_p\) was correct, but its paired shared-host stage time
was essentially level with the original routine. Keeping the original
extension-field interpolation while moving only the resultant to
\(\mathbb F_p\) was slower in this run. These results direct the next
optimization experiment toward the fiber-lift and interpolation phases.

| Same-input method | Exact cuts matched | Sum of 32 stage intervals | Median interval per cut | Paired wins versus original |
| --- | ---: | ---: | ---: | ---: |
| Original extension-field routine | 32/32 | 24.322 s | 671.696 ms | Reference |
| Prime-field coefficient rows and resultant | 32/32 | 24.467 s | 741.320 ms | 17/32 |
| Extension-field interpolation, prime-field resultant | 32/32 | 26.532 s | 781.473 ms | 11/32 |

The input is the first 16 rows of the [frozen two-pivot reduction](../../type_i_two_pivot_reduced64.json),
with both quotient point `Q` and its unique odd-order half `H=Q/2` checked
for each row. The [producer](../../line_quotient_cubic_fast_probe.py)
constructs the same inverse-cover fibers and plane-quartic coordinates as
the existing cubic-cut constructor. The coefficient-row variant evaluates
the ten cubic monomials at two fiber points and three rational anchors,
expands each value in the six-dimensional \(\mathbb F_p\)-basis, and
solves the resulting prime-field linear system. The six geometric fiber
intersections are recovered as the two three-step Frobenius orbits in the
known-intersection polynomial. It then computes the quartic/cubic
resultant and linear subresultant over \(\mathbb F_p\). The hybrid keeps
the original nine-point extension-field cubic interpolation and moves
only those elimination steps to \(\mathbb F_p\). Both methods assert
equality of all cubic coefficients, residual `x` polynomials, and residual
`y` polynomials against the frozen GP-verified receipt for every cut.

Each method is timed only while constructing one residual cut after
the common Sage process and model setup. The three method orders rotate
by class index, and each class supplies both `Q` and `H` cuts. The
[producer receipt](../../type_i_fast_cubic_probe.json) retains all 96
individual intervals and the exact input/source hashes; the
[independent receipt audit](../../type_i_fast_cubic_probe_verification.json)
recomputes stage totals, source binding, runtime binding, and progression.
The checked [Sage runtime receipt](../../type_i_fast_cubic_runtime_info.json)
was saved before the experiment. These are unisolated implementation-stage
measurements, so the near-equal totals do not support a device-level
performance conclusion. The earlier single-class probes and their receipts
are preserved as attempts; they used an unbalanced method order.

The coefficient-row equivalence has a short Lean proof:
`P192ExtensionCover.cubic_evaluation_zero_iff_coordinate_rows_zero` in
[CubicCoefficientRows.lean](../../lean/CubicCoefficientRows.lean) says
that an extension-field evaluation vanishes exactly when all coordinates
vanish, under an injective linear coordinate map. The concrete six-term
P-192 basis, smooth quartic, fiber points, and polynomial intersections
are checked by the Sage producer and the previously frozen native GP
replay. The Lean theorem is a conditional algebraic bridge, not a
formalization of those concrete geometric instances.

For the dangerous-curve screen, this experiment leaves the exact
three-sum support and known-scalar query law from the preceding run
unchanged: the tested 48-class base has 4,926 admitted targets and zero
rank-bearing triple collisions. The full 20,014 signed fan columns are
available as exact quotient coordinates. Materializing all their reduced
divisors with this cubic-cut path would expand representation coverage,
while the existing fixed-base support bound still governs this
three-sum policy. Four of the 20,014 distinct two-fan columns coincide
with selected classes. Granting all 9,925 selected classes and every
remaining fan column as usable base points gives `B=29,935` and at
most `C(B+2,3)=4,471,261,391,320` three-sum targets. Against the
certified accepted scalar set, the ideal uniform accepted-query hit
probability remains below `2^-147`. The additional Lean declarations
`full_two_pivot_base_three_sum_capacity` and
`full_two_pivot_base_hit_bound_certificate` check those integers; their
finite-field distinctness premise is supplied by the two-fan certificate.
The [capacity receipt](../../type_i_two_fan_capacity.json) reconstructs
the full exact quotient-point union and reads the accepted-set lower
bound from the frozen density receipt. A separate
[sorted-list verifier](../../type_i_two_fan_capacity_verification.json)
reproduces both two-point overlaps, the 29,935-class union, and the
capacity without using the producer's sets. The pinned
[`PivotGuided` build](../../lean/pivot-guided-full-fan-build-20261016.log)
checks the integer implications with an axiom audit; the
[`CubicCoefficientRows` build](../../lean/cubic-coefficient-rows-build-20261016.log)
checks the conditional coordinate-row theorem.
The next profiling gate is an exclusive breakdown of
fiber lifts, cubic interpolation, known-intersection construction, and
resultant/subresultant costs on a frozen wider sample.

Reproduce from the `crypto` worktree root:

```sh
/Volumes/SSD990/cryptanalysis/sage --runtime-info
/Volumes/SSD990/cryptanalysis/sage -python research/p192_extension_cover_20261012/line_quotient_cubic_fast_probe.py --count 16
python3 -B research/p192_extension_cover_20261012/verify_fast_cubic_probe.py
cd research/p192_extension_cover_20261012/lean && lake build CubicCoefficientRows
```

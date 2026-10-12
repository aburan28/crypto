# Stage 191: selected F4 single-core and direct-MITM boundary

## Result

The selected Stage 190 F4 source remains far slower than same-binary direct
meet-in-the-middle decomposition on the identical opened public
`n=59, ell=9, m=3` true-negative.

The fixed fresh-process schedule was direct MITM, strict single-core F4, direct
MITM. The direct reference is the arithmetic mean of both bracketing controls;
neither the faster control nor an inherited timing is selected.

| arm | wall seconds | total core-seconds | single-core seconds | peak RSS |
|:---|---:|---:|---:|---:|
| direct MITM before | 0.746381 | 0.308969 | 0.308969 | 41,222,144 B |
| direct MITM after | 0.324357 | 0.306830 | 0.306830 | 41,517,056 B |
| **direct arithmetic-mean reference** | **0.535369** | **0.307900** | **0.307900** | **41,369,600 B** |
| **selected F4, strict single core** | **109.175286** | **109.121473** | **109.121473** | **480,526,336 B** |

The resulting F4/direct ratios are:

- wall: `203.925343`;
- total core-seconds: `354.406139`;
- single-core seconds: `354.406139`; and
- peak RSS: `11.615446`.

Direct wall has visible host-noise variation, while its CPU values are stable.
The frozen arithmetic mean keeps both controls in the reference. The CPU ratio
is the more stable boundary diagnostic.

## Correctness and work

Both direct controls authenticate the same source, materialize 483 rational
factor-base points, store 116,635 pair entries, charge 117,369 group additions,
and return exhaustive `UNSAT` with no witness.

The F4 arm uses `RAYON_NUM_THREADS=1`, `PQ_F4_X1_BATCH=1`, all external thread
caps set to one, dense exact pair selection, full M4RI disabled, and no Stage
190 scheduling override. It reports `single_thread_requested=true`, 242
selected build-serial calls, all 512 masks visited, 270 non-rational masks
skipped, all 242 systems completed, 14,278 equations, 14,515,915 terms, zero
roots, and exhaustive `UNSAT`.

The F4 work counters remain 319,313,687,585 row-equivalent XORs and
147,794,583,858 actually performed/table-construction XORs. Both methods use
the algebraic `span_F2(1,z,...,z^8)` base without target-subgroup enumeration
or known discrete-log labels. Conflicts are `null` because neither method
exposes a SAT-conflict counter.

Final native replay passes `24/24` with result SHA-256
`7708914c8e1ac5a0dfec6bd5b4bdb3fed32669aa17b1e189769c10dad8ba781d`.

## Accounting and claim boundary

Stage 191 contributes a measured lower bound of 8 components,
`237.736168` wall-seconds, `655.247344` total core-seconds, and
`5,625,905,152` bytes peak RSS. The cumulative measured campaign lower bound is
572 components, `23,714.619162` wall-seconds, `60,598.525724` core-seconds,
and `6,310,576,128` bytes maximum RSS.

Complete campaign cost remains `null`; development compilation before the
exact native build and final writes after the charged composition control are
not all outer-metered.

This is a decomposition-stage comparison. It does not price natural relation
yield, relation collection, relation linear algebra, target descent, log
recovery, or the full automorphism-aware rho reference. Gate 6 remains false,
gate 7 remains open, and this is not Koblitz index-calculus SOTA.

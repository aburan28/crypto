# Exact layout reuse and fused packing for inherited-F4 roots

## Result

Inherited-F4 roots now reuse a column layout across public targets only when
the cached support is exact. On a hit, shifted polynomial rows are formed,
cancelled, and packed directly into the nested dense matrix. A missing product
or an unobserved cached column rejects the hit and rebuilds the support.

The automatic policy applies to quadratic generator systems. The cubic frozen
cell produced seven hits and eight support misses in the exploratory panel, so
it stays on the historical uncached builder. `KIC_F4_INHERIT_LAYOUT_CACHE=1`
retains forced experiments, `=0` is the uncached control, and
`KIC_F4_DISABLE_INHERIT_FUSED_PACK=1` keeps materialized sparse rows on exact
layout hits.

All nine nonignored F4/RREF tests pass. They include exact cached-versus-
uncached root matrices, deliberate changed-support fallback, nested and flat
fused packing at degrees 2 through 4, missing/extra-support rejection, and the
existing RREF reference checks.

## Root packing

Each public-synthetic Semaev root has 100 alternating materialized/fused
samples after the exact layout is fixed:

| fixture | shape | materialized median | fused median | speedup |
|:--|--:|--:|--:|--:|
| binary n=13, ell=12, m=2, d=3 | 325 x 1,885 | 0.198 ms | 0.175 ms | **1.129x** |
| binary n=19, ell=18, m=2, d=3 | 703 x 6,175 | 1.172 ms | 1.100 ms | **1.065x** |
| binary n=23, ell=11, m=2, d=3 | 529 x 1,464 | 0.309 ms | 0.275 ms | **1.122x** |
| binary n=31, ell=16, m=2, d=3 | 1,023 x 4,369 | 0.720 ms | 0.645 ms | **1.116x** |

Every matrix is byte-identical. The aggregate paired-ratio median is
**1.115x**.

## Whole-stage panel

Ten fresh-process repetitions per policy compare automatic cached+fused roots,
automatic cached+materialized roots, and fully uncached roots. All verdict
digests, decomposition counts, solver trees, matrix shapes, oversize counters,
and word-operation counts match exactly.

| panel | fused vs materialized stage | fused vs uncached stage | fused vs uncached child CPU | fused vs uncached build |
|:--|--:|--:|--:|--:|
| frozen | **1.000x** | **1.069x** | **1.065x** | **1.127x** |
| holdout | **1.009x** | **1.173x** | **1.024x** | **1.268x** |

Every quadratic rung improves against uncached construction, from **1.028x**
through **1.242x** on the frozen ladder and **1.061x** through **1.242x** on
holdout. The cubic rung records zero layout hits and misses under the automatic
policy and therefore preserves its old path.

## Boundary

This is public-synthetic decomposition-stage engineering evidence. It changes
root construction overhead only. It does not change the factor base, Macaulay
row space, degree, solver tree, relation yield, asymptotic exponent, or the
candidate-tuple floor. It is not an end-to-end index-calculus crossover, a
sub-Pollard-rho result, an external-target capability, or a key-recovery result.

# Reusable M4RI arenas for inherited matrix-F4

## Result

The block-4 Method of Four Russians kernel now retains its pivot-column list,
pivot row buffers, and combination table in one thread-local arena. The prior
implementation allocated all three again for every pivot block. The original
implementation remains selectable with `KIC_F4_M4RI_ALLOCATING=1`, so every
measurement below is a same-binary comparison.

The exact gate passed 96 combinations of matrix shape, density, seed, block
width (2 through 5), and reduction mode. Arena and allocating implementations
returned byte-identical matrices, identical ranks, and identical counted word
operations in every case.

On 600 alternating synthetic samples, every cell improved. Median speedups
ranged from **1.095x to 1.458x**, with an aggregate paired-ratio median of
**1.201x**. Both plain echelon form, used by wide inherited roots, and full
RREF, used by narrower roots and ordinary matrix-F4, passed.

The real public-synthetic Semaev root matrices are the decisive control:

| fixture | shape | root form | allocating median | arena median | speedup |
|:--|--:|:--|--:|--:|--:|
| binary n=13, ell=12, m=2, d=3 | 325 x 1,885 | echelon | 0.240 ms | 0.216 ms | **1.114x** |
| binary n=19, ell=18, m=2, d=3 | 703 x 6,175 | echelon | 1.956 ms | 1.647 ms | **1.188x** |
| binary n=23, ell=11, m=2, d=3 | 529 x 1,464 | full RREF | 0.590 ms | 0.537 ms | **1.099x** |
| binary n=31, ell=16, m=2, d=3 | 1,023 x 4,369 | echelon | 2.546 ms | 2.158 ms | **1.180x** |

Each row contains 100 alternating samples after warmup. The aggregate paired
ratio median is **1.131x**. Ranks, complete output matrices, and operation
counts match on every sample.

## Whole-stage control

The current default inherited-F4 engine builds a Macaulay matrix only at each
target's root and specializes the retained basis below it. Root elimination is
therefore a small part of the already-optimized stage. Across 13 fresh-process
frozen repetitions per policy, child CPU favored the arena by 1.033x while
stage wall favored the allocating control by 1.050x. Across 15 holdout
repetitions, both were within about 1.4% and slightly favored the control.
These conflicting host-sensitive metrics are retained as **indeterminate**,
not promoted as a whole-stage wall-time win.

Every frozen and holdout pair has identical verdict digests, decomposed target
counts, reductions, refutations, propagations, splits, matrix calls, rows,
columns, oversize counters, criterion costs, specialization costs, and total
word operations. The optimization changes allocation overhead only.

## Boundary

This is public-synthetic decomposition-stage engineering evidence. It does not
change the Macaulay row space, splitting tree, relation yield, factor base,
degree of regularity, asymptotic exponent, or the candidate-tuple floor. It is
not an end-to-end index-calculus crossover, a sub-Pollard-rho result, an
external-target capability, or a key-recovery result.

# Active-variable and exact-layout F4 construction

This iteration removes construction work that the specialised Boolean system
cannot use. Once a variable disappears from every generator, multiplying by it
cannot create a consequence without that variable. The solver therefore uses
only currently occurring variables for Macaulay shifts on systems with at
least 24 variables.

Multiplier schedules are cached per thread by active mask and degree. Column
layouts are cached by active mask and Macaulay degree, but every hit is accepted
only after every row monomial maps into the cached layout and every cached
column is observed. A mismatch rebuilds and replaces the layout. Both caches
are bounded and have same-binary disable controls.

## Final same-binary pair

| Measurement | Merged-policy control | Selected | Ratio |
|---|---:|---:|---:|
| Complete process wall | 150.285 s | 114.041 s | 1.318x |
| Child user CPU | 149.767 s | 113.192 s | 1.323x |
| x-chained found median | 27.502 s | 12.203 s | 2.254x |
| Symmetrised refuted median | 15.331 s | 7.441 s | 2.060x |
| Peak RSS | 467,206,144 B | 485,965,824 B | 1.040x selected/control |

Both policies used the same executable and two public targets. Both x targets
remain found, both symmetrised targets remain refuted, effort remains
7,056/4,095, built degree remains 3, FFD remains 3/4, every gate passes, and
there are no inconclusive outcomes.

Relative to the pre-linear-tail legacy retained by the predecessor package,
the cumulative selected process is about 1.9x faster; the x and sym arm medians
are more than four times faster.

## Mechanism

The sym selected profile takes 7.323 s versus 15.238 s control. Rows fall from
8,125,472 to 5,331,938, accumulated columns from 13,561,904 to 6,524,366, and
word XORs from 4,509,530,169 to 2,415,432,155. It records 8,178 exact layout
hits and 13 misses. Disabling layout reuse raises wall to 8.707 s; disabling
the small multiplier-schedule cache raises it to 7.391 s.

The x first-root profile takes 0.895 s versus 2.002 s control, with 1,013 exact
layout hits and 15 misses. Its XOR count falls from 584,547,728 to 265,653,883.

All frozen and holdout verdict digests match. The three above-gate ladder rows
improve 2.707x, 3.023x, and 1.610x. Below-gate rows keep identical operation
counts and run within normal timing noise of the control.

## Rejected mechanisms

- A dense-row allocation arena is about 11% slower because retained capacities
  harm reduction locality.
- Verified pivot-trace replay with a fallback produces no measurable reduction
  gain; validating skipped ranges costs about what their scans save.
- Hash-based column collection builds more slowly than sort-and-deduplicate.

Classification: `N31_BOOLEAN_F4_ACTIVE_LAYOUT_PAIRED_T2_PASS`.

This is bounded public-synthetic decomposition-construction evidence. It does
not measure relation yield, a complete index-calculus run, asymptotic
complexity, or a crossover against Pollard rho, and it is not a key-recovery
result.

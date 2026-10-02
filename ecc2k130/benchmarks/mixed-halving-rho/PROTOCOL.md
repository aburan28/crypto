# Point-dependent mixed-halving rho: frozen finite screen

Baseline source commit: `1c2e8e1898ac89fb8b3e8880761e4a2eda196df2`.

This is an exact host-only functional-graph screen. It neither runs the
ECC2K-130 search nor claims GPU throughput. The purpose is to decide whether a
point-dependent mixture deserves a CUDA prototype after the globally scheduled
halving design failed its collision-normalized model.

## Candidate map

Let `d = 2^-1 mod ell`. A Frobenius- and negation-invariant selector word is
computed from the current point. Its low four bits select one of these halving
fractions:

```
0, 8/16, 12/16, 14/16, 15/16, 16/16.
```

If the low nibble is below the threshold, the map takes the unique odd-subgroup
half `[d]R`. Otherwise it uses bits 4..6 to select one of eight existing
equivariant table addends and takes `R + A_h(R)`. Unlike a global schedule,
the choice is a function of `R`; there is no phase coordinate and no family of
correlated intermediate states.

Two selector words are tested:

- `canonical_hash`, an expensive ideal-balance control formed from the full
  Frobenius-orbit canonical x coordinate;
- `necklace7`, the seven inexpensive rotation invariants from the stateless-v4
  screen, mixed with normal-basis x weight before the bit split.

The 0/16 row is the additive structural-cycle control. The 16/16 row is the
permutation control. Neither can pass. At the measured primitive rates
`R_half=28.953 B/s` and `R_add=20.1343255 B/s`, 12/16, 14/16 and 15/16 have
ideal free-dispatch raw rates above 26 B/s. Those rates omit representation,
dispatch, distinguished-point and guard costs and are not results.

## Exact experiment

For the generated degree-23 prime subgroup (`ell=2,095,853`), enumerate all
nonidentity points and all 45,562 classes under Frobenius and negation. For 16
frozen target/table seeds and every selector/threshold pair:

1. Verify the selector, phase, sign and next-class covariance exhaustively.
2. Build the complete point and quotient functional graphs.
3. Record indegree collisions, exact cycles, the basin ending in cycles of
   length at most eight, and each state's tail-plus-cycle first-repeat length.
4. Compare with 16 deterministic random functional graphs of the same quotient
   size.

The first-repeat statistic is a finite diagnostic. It is not substituted for a
challenge collision constant.

## Admission rule

A CUDA prototype is admitted only if a row with an ideal raw rate above
26 B/s satisfies all of these pooled over the 16 seeds:

- zero covariance failures;
- actual halving fraction within 0.02 of its threshold fraction;
- nonzero quotient indegree collisions on every seed;
- mean quotient two-cycle count no greater than twice the random-control mean
  plus one;
- mean quotient short-cycle-basin fraction no greater than twice the
  random-control mean plus 0.02;
- median across seeds of `mean_first_repeat / sqrt(N)` between 0.75 and 1.50
  times the random-control median;
- the eight addition branches have pooled `sum(p_h^2)` no more than 1.10 times
  uniform.

Passing admits only a representation and dispatch cost study. A CUDA walk must
then preserve coefficient replay, point canonicalization, finite guards and
distinguished-point reporting, and must beat 26 B/s end to end before the goal
is met.


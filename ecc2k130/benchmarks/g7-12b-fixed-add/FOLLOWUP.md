# Preregistered follow-up: Frobenius-covariant fixed add

Status: closed at the preregistered CPU pilot before CUDA implementation. With
a 1,024-step restart, mean collision work was 1.183749x selected, above the
1.10 continuation threshold. The default no-restart mode also has a witnessed
four-state DP-free cycle.

## Hypothesis

The first fixed-add candidate is fast because it removes two normal-to-polynomial
conversions, one polynomial-to-normal conversion and both runtime Frobenius
routes per point. It loses collision efficiency because a fixed branch addend
does not commute with Frobenius.

For a nonexceptional normal-basis coordinate `x`, number the coordinate positions
along the Frobenius cycle and define

```
phase(x) = sum(position of each set bit) / weight(x) mod 131.
```

Frobenius advances every set-bit position by one, so `phase(sigma(x)) =
phase(x)+1`. Keep the Hamming-weight branch `j=3..10`, but add
`sigma^(phase(x)+j)(B)`. Then `T(sigma(P))=sigma(T(P))`; the walk again descends
to the Frobenius quotient while its point formula still consumes polynomial
fixed-point coordinates. For weight zero or 131, use the selected self-Frobenius
step because phase is undefined; this branch is also covariant.

Replay tracks both scalar coefficients directly. An additive step increments the
B coefficient by `s^(phase+j)`; an exceptional self-Frobenius step multiplies
both coefficients by `1+s^j`.

## Gates and decision

Before CUDA work, independently verify the phase identity for every
nonexceptional GF(2^23) coordinate and run the same planted-collision harness as
the selected/fixed-add study. Use identical seeds, DP rule, partition and orbit
key. First run 100 trials per mode. Continue to 2,000 only if every scalar is
recovered and mean complete collision iterations are no more than 1.10x
selected. Close the candidate if the 2,000-trial upper evidence remains above
1.10x. Otherwise implement it in a fresh isolated source copy.

The CUDA candidate will store all 131 polynomial orbit addends in shared memory;
divergent constant-memory indexing is disallowed. It must prove exact phase,
addend, full-state, DP, replay and collision behavior against an independent CPU
reference before timing. The initial screen is three interleaved equal-work
selected/candidate pairs on the same g7.2xlarge / RTX PRO 4500 at 165 W. Report
both raw and collision-adjusted paired log-ratio intervals. Promotion still
requires both intervals wholly above one, followed by a fresh confirmation.
Generic work is `sqrt(n/262)` times the measured collision-work ratio. Full-DLP
S remains null until measured end to end.

## Harness correction before completed pilot

The first 100-trial invocation checked the scalar-coordinate invariant after
every distinguished walk and was interrupted after four minutes with no result;
the field scalar multiplications dominated the collision measurement. Preserve
the empty partial output. Check both invariants only when an orbit-key collision
occurs, which still gates every attempted scalar solve without perturbing paths
or work counts, then restart the pilot from trial zero.

The restarted pilot reproduced both control rows, then exposed avoidable harness
cost: it recomputed `sigma^phase(B)` inside every CPU step. Interrupt before the
covariant result, preserve the partial output, precompute the 23-point orbit just
as the CUDA proposal precomputes 131 addends, and restart all rows. This changes
neither paths nor reported iteration counts.

The precomputed-orbit restart again reproduced both controls but the covariant
row did not finish within four minutes. Preserve that partial output. Add a
mode-only diagnostic argument and trial progress, with no path or result change;
run one covariant trial before committing more CPU time.

The isolated 100-row run reached the final ten trials and then exceeded the
one-trial diagnostic by several minutes. Interrupt at 7m42s and preserve its
partial time/progress logs. Add a trial-offset diagnostic argument and report
every trial, without changing the walk, to isolate the rare case before using an
aggregate.

The offset diagnostic identifies trial 94 as the stall. Covariant additive walks
preserve Q coefficient one and frequently coalesce at the exact same point;
`c=0, epsilon=+1` makes the collision denominator zero and supplies no discrete
log equation. Count and skip these known-degenerate exact duplicates before the
expensive scalar-coordinate invariant check. The prior harness already skipped
them after computing the same zero denominator, so this changes cost only.

Trial 94 still stalled after exact-duplicate rejection, showing it is trapped
inside one walk rather than processing collisions. Interrupt and preserve the
partial output. Add a 100,000-step diagnostic cap that records the trial, walk
and deterministic seed, then continues with a new seed. This is diagnostic
restart accounting; any hit is a candidate defect because the production
default currently disables `--max-iters`.

The exact-state witness proves the offending seed begins in a four-state cycle;
its minimum x weight is 12 while the DP threshold is 10. No exceptional phase
branch is involved. Evaluate the existing restart mechanism at 1,024 steps and
charge every overdue step to collision work. This is the smallest production-
compatible repair considered; CUDA work still requires the repaired 100-trial
mean to stay within 1.10x selected with zero bad recoveries.

Correction: exact-point collisions are only known-degenerate when their scalar
Q coefficients are also equal. The additive candidate has `b=1`, but selected
self-Frobenius walks do not. Restrict the fast rejection to zero `b_B-b_A` so
the control retains its informative exact collisions. This does not change the
completed covariant row, whose exceptional-step count is zero.


## Result

The final matched 100-trial control solved every scalar with 166.640 mean
complete iterations. The covariant candidate with a 1,024-step restart solved
every scalar with 197.260 mean complete iterations, one overdue restart and
three zero-denominator exact duplicates. Its mean-work ratio is **1.183749x**,
so the preregistered 1.10 gate closes it without a 2,000-trial study or CUDA
implementation.

The deterministic seed `5655958182028705792` begins in an exact four-state
cycle. Its minimum x-coordinate weight is 12, above the GF(2^23) DP threshold
10. Production's default `maxIters=0` would never escape it. With restart 1024,
generic work is `sqrt(n/262) * 1.183749`; full-DLP S remains null.

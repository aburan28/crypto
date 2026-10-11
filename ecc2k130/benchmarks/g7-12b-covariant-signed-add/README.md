# Signed covariant fixed-add rho toward 12 B/s on G7

Status: preregistered before implementation, compilation, collision study or GPU timing.

## Hypothesis

The non-equivariant fixed-add kernel reached 7.601382 billion raw updates/s but
lost collision efficiency. Its phase-covariant successor restored Frobenius
covariance but not negation covariance, averaged 1.183749x selected collision
work with restart1024, and exposed a four-state DP-free cycle.

For nonexceptional normal-basis x, retain the centroid phase

```
phase(x) = sum(Frobenius-cycle positions of set bits) / weight(x) mod m.
```

Normalize each set-bit position by subtracting phase, choose the least normalized
position, and read y at its corresponding set x coordinate. Call that bit
`sign(P)`. Frobenius shifts x, y and phase together, so sign is invariant.
Negation changes y to y+x, so the chosen y bit flips. Let `A=B+Q`; branch j adds
`sigma^(phase+j)(A)` when sign is zero and its point negative when sign is one.
The resulting transition satisfies both

```
T(sigma(P)) = sigma(T(P))
T(-P)       = -T(P).
```

Using A=B+Q changes both replay coefficients by plus or minus `s^(phase+j)`.
Exact point collisions therefore need not retain equal Q coefficients. Weight
zero or m uses the selected self-Frobenius transition. Restart at 1,024 steps
and charge all overdue work.

## Frozen gates

First exhaustively prove phase and selected-bit covariance for every
nonexceptional GF(2^23) x coordinate. Prove sign covariance and negation flipping
on at least 131,072 valid seeded points. Then run 100 matched planted GF(2^23)
DLPs for selected and candidate with identical seeds, DP rule, partition, orbit
key and restart1024. Every scalar and every checked full-state coefficient must
recover. Continue to 2,000 trials only if candidate mean complete collision work
is at most 1.10x selected and no unexplained cycle, infinity or invariant failure
occurs.

Only a passing collision study permits an isolated CUDA implementation. Store
all 131 polynomial A-orbit points in shared memory, compute phase/sign from the
normal x/y already needed by the selector, and keep the production feature off
by default. Exact state, DP, replay, resume, restart and sanitizer gates precede
three interleaved equal-work selected/candidate G7 pairs at 165 W. Report raw and
collision-adjusted paired log-ratio intervals. Promotion requires both intervals
wholly above one and a fresh confirmation. The measured target remains 12
billion complete iterations/s. Generic work is `sqrt(n/262)` times the measured
collision-work ratio; full-DLP S stays null until an end-to-end recovery exists.

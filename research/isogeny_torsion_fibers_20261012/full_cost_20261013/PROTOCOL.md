# Frozen growing-size, complete native-cost comparison

Date: 2026-10-13. Conductor T-2486. This protocol is committed before
candidate selection or solver execution. It follows the verified F617
[non-isomorphic pair](../nonisomorphic_pair_20261013/PROTOCOL.md) and retains
that run as historical evidence.

## Question, boundary, and fixed inputs

Measure whether transporting each public source target through the rational
degree-three isogeny improves a complete cold relation solve. Select one pair
in each inclusive prime-subgroup interval `[101,400]`, `[401,1600]`,
`[1601,6400]`, `[6401,25600]`: the smallest prime `r` in the interval for
which `p=6r-1` is prime. The first must be `(r,p)=(103,617)`.
For each, use `E/F_p: y²=x³+1`, `E'/F_p: y²=x³-27`, subgroup order `r`,
cofactor six, and

`phi(O)=phi(0,±1)=O; phi(x,y)=((x³+4)/x², y(x³-8)/x³)` for `x≠0`.

Both models have order `p+1=6r`. Verify by enumeration. Enumerate sixth
powers to reject an `F_p` scaling isomorphism. Check every source point's
image is on the destination, the rational kernel is exactly three points,
and the image fibers each have size three. Check the homomorphism exhaustively
at F617 and on a fixed 256-pair sample plus every subgroup multiple at the
larger sizes. Any structural failure stops solver execution and is retained.

The fixed generic two-class walk-step floor is `sqrt(pi/4)` in
`S=(group additions+doublings)/sqrt(r)`. The measured exact-target reference
is source negation rho. The primary complete-cost unit is **cold elapsed
nanoseconds per target and arm**, timed from before method-specific setup
through verified recovery, including the isogeny evaluations. This direct
unit charges roots, table operations, row arithmetic, hashes and allocation
without assuming a conversion to group additions. Report group-operation `S`
only as a diagnostic, alongside native counters and peak memory. Full elapsed
time is admitted for a speed claim only if an isolated L2-or-better session,
an A/A control, and an independent host replay are obtained. A shared-host
L0 run is a complete-cost observation, with its isolation limit stated.

## Target and arm contract

For each pair choose the first affine source point in `(x,y)` order whose
`[6]` projection is nonzero as generator `G`; set `G'=phi(G)` and check both
orders. Generate eight distinct nonzero public subgroup targets by the F617
SplitMix64 lift-and-project rule, seed `0x2026101300010000+p`. Freeze exact
coordinates in raw output; no target scalar enters a solver.

Each target has five interleaved measured rounds and one warmup round. Every
round has source relations, transported destination relations, matched source
negation rho, and a duplicate rho A/A arm. The duplicate uses the same seed
and algorithm; its count/output identity and timing spread measure noise.
Rotate arm order by `(round+target_index) mod 4`. Each relation arm builds
its factor base and pair table anew; both use twelve projected points,
the same trial streams and `100000` collection/descent caps. Rho has a
`1000000` step cap. Seeds are the F617 seeds XORed with `round<<32`;
all arms on the same target/round share their stream keys. An exhausted
attempt is recorded and charged. No best-of selection or retry is allowed
without a separately frozen version.

The destination arm times both evaluations `phi(G)` and `phi(Q)`. All
relation arms time factor-base lifts, projection, pair-table construction,
collection, linear algebra, descent, and scalar verification. Rho times
class setup, walk, collision recovery and verification. Every returned
scalar is checked on the original source; destination rows are checked on
both curves. Target generation, structural checks and candidate search are
reported as construction/preflight separately because a public ECDLP input
already supplies the curve, generator and target. Report a second cold-start
cost including construction; never amortize setup across targets in the
primary row. Record process RSS, phase time and counters, including maps,
roots, lookups, modular work, and failures. Compare sums over all attempts;
bootstrap paired target intervals with a fixed seed and retain the raw rows.

## Validation and limits

Build and run native Rust only. Require all 4×8×5×4 measured rows to have
verified equal answers; retain any partial or failed raw version. Run a
byte-identical native deterministic replay of the answer/count records and
an independent Rust arithmetic checker. Timings are never expected to
replay byte-for-byte. Run the required release library and example tests.
Record source/binary hashes, exact commands, host, isolation, elapsed phase
totals, memory, native counters, paired ratios and intervals, a visual and
PDF. Copy performance visuals and PDFs to `benchmarks/`.

If the local host cannot earn L2 or no independent runner is available,
retain and report the complete-cost L0 observations and keep the admitted
speed comparison open. A bounded failure or timing difference alone does not
establish a general scaling result or an advantage on another curve.

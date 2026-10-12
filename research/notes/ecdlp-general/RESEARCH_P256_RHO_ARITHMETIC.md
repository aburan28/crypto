# P-256 Pollard-rho arithmetic experiment plan

**Status:** pre-registered; no performance result is claimed by this note.

**Question.** Can P-256's generalized-Mersenne field representation reduce the
*end-to-end* cost of a Pollard-rho iteration when arithmetic, normalization,
partitioning, table walks, batching, and the target architecture are designed
together?

The prime is

```
p = 2^256 - 2^224 + 2^192 + 2^96 - 1.
```

This is an engineering thread, not a claim of a sub-generic attack.  The
generic-group reference remains Pollard rho on the same subgroup.  All
candidate variants must solve the same planted DLP instances and are compared
to the same unmodified rho implementation.

## Boundary and accounting

For a subgroup of order `n`, define

```
S = counted curve/field-operation cost / sqrt(n)
C = observed collision steps / sqrt(pi*n/2)
R_eff = raw_steps_per_second / C
```

The repository's primary metric remains counted end-to-end cost `S`; wall
clock and `R_eff` are engineering diagnostics.  A candidate is not called a
speedup unless the complete solve, including table construction, restarts,
cycle escape, distinguished-point handling, collision recovery, and
verification, is cheaper than the baseline under the same accounting.

**Reference:** the existing r-adding rho implementation with its current
cycle handling and distinguished-point machinery, adapted to P-256 without
the secp256k1-only order-3 automorphism.  P-256 may use fold 1 or the negation
map (fold 2); this experiment does not assume GLV/Aut(E)=Z/6.

## Frozen inputs

Before timing, freeze:

- curve: NIST P-256, plus toy short-Weierstrass curves for complete DLP solves;
- toy subgroup sizes: 32, 36, 40, and 44 bits where practical;
- seeds: `0..63` for toy complete solves; field microbench seed `0x50323536`;
- r-adding table sizes: `R = 8,16,32,64,128,256`;
- normalization intervals: `K = 1,2,4,8,16,32`;
- simultaneous-inversion widths: `B = 8,16,32,64,128,256,512,1024,2048`;
- DP densities: `2^-12, 2^-16, 2^-20, 2^-24, 2^-28`;
- identical table coefficients, starting points, secrets, and DP rules across
  matched variants.

Every generated fixture must record the generator version/hash and seed.

## E1: field-backend shootout

Implement interchangeable P-256 field backends:

1. Montgomery CIOS;
2. P-256-specialized Montgomery;
3. eager NIST/Solinas reduction;
4. fused multiply + P-256 reduction.

Benchmark random multiplication, squaring, dependent multiply chains, and the
actual field-operation sequence emitted by point addition.

Record cycles/op, multiplies, adds/subtracts, carry-chain instructions,
registers/thread, spills, occupancy, and throughput where available.

**Hypothesis:** fused reduction lowers instruction count, but the rho winner
may differ from the isolated field-multiply winner because of register
pressure and carry dependencies.

**Success:** candidate lowers complete rho cost on at least two toy sizes and
does not regress correctness.

**Stop:** reject a backend if it loses >=10% end-to-end after architecture
tuning even if its isolated multiplication benchmark wins.

## E2: fused Solinas multiply/reduce

Instrument the product/reduction boundary.  Compare:

- materialize 512-bit product then reduce;
- interleave high-limb folding with multiplication;
- architecture-specific schedules that consume high limbs as soon as their
  final contribution is known.

Differentially test every output against a big-integer oracle on edge cases
`0,1,p-1`, limb carry/borrow adversaries, and at least 100,000 seeded random
pairs.

**Required invariant:** exact equality modulo P-256 for every test vector.

Report both isolated latency and the change in the complete point-add/rho
kernel.  A stage-only improvement is labelled a stage diagnostic.

## E3/E4: weak reduction and normalization-depth sweep

Represent field values with an explicit redundancy bound.  Start with
`[0,p)`, `[0,2p)`, `[0,4p)`, and `[0,8p)`; larger ranges are admitted
only after an overflow proof/check.

For each field primitive, record the maximum observed limb magnitude and
cross-check against arbitrary-precision arithmetic.  Then sweep forced
normalization every `K = 1,2,4,8,16,32` rho iterations.

Measure:

- canonicalizations avoided;
- carry/borrow operations;
- maximum redundant range reached;
- raw steps/s;
- collision constant `C`;
- complete solve cost.

**Hypothesis:** an intermediate `K` wins; "never reduce" is not assumed.

**Stop:** any unproved overflow range or oracle mismatch immediately rejects
that representation.

## E5/E6/E7: partition-function cost, quality, and invariance

Baseline:

```
h0(P) = x_low & (R-1)
```

Candidates use only already-live limbs, for example:

```
h1(P) = (x0 ^ x1) & (R-1)
h2(P) = mix(x0,x1) & (R-1)
```

No candidate may require a canonical reduction merely to choose the next
table entry.

Three gates are mandatory:

1. **cost:** cycles/instructions per partition;
2. **quality:** histogram bias, serial correlation, short-cycle rate, escape
   rate, and complete-solve collision constant;
3. **representation invariance:** equivalent representations of the same
   affine point must select the same partition whenever the walk state admits
   multiple encodings.

For Jacobian coordinates, test random nonzero `lambda` under
`(X,Y,Z) -> (lambda^2 X, lambda^3 Y, lambda Z)`.  A raw-projective hash that
changes under this transform is rejected as a group-element partition.

**Success:** lower partition cost with collision constant within the
pre-registered confidence interval of baseline and no increase in pathological
cycle rate.

## E8: r-adding table sweep

Sweep `R = 8..256` powers of two using identical seeded table generation.
For each R, run the complete toy solve panel and record:

- raw steps/s;
- mean/median collision steps;
- `C = steps/sqrt(pi*n/2)` (or the fold-2 denominator when enabled);
- short-cycle detections by length;
- escape steps;
- abort/reseed rate;
- table-cache behavior.

Do not select R on raw throughput alone.  Compare complete expected solve cost
and preserve all seeds, including slow ones.

## E9/E10: Jacobian versus batched affine

Compare three matched walkers:

A. Jacobian/mixed-add;
B. affine with one inversion per walker-step;
C. affine with Montgomery simultaneous inversion over B walkers.

Sweep `B = 8..2048`.  The table, partition, starts, secrets, and DP policy
must be identical.

Measure inversion-equivalent field cost, temporary storage, registers,
occupancy, raw steps/s, `C`, and complete solve cost.

**Hypothesis:** enough independent walkers can make batched affine cheaper than
a conventional projective hot loop.

**Critical correctness rule:** the partition remains a deterministic function
of the group element.  Projective scaling must not alter the walk.

## E11: distinguished-point batching

Sweep DP densities `2^-12 .. 2^-28`.  Separately measure DP test,
normalization, emission, host transfer, duplicate handling, replay/collision
recovery, and restart cost.

Compare normalization policies:

- every step;
- every batch;
- only at a representation boundary needed by partition/DP handling.

The DP test must not silently force a full affine/canonical conversion every
iteration unless that cost is explicitly counted.

## E12: CPU matrix

At minimum compare 4x64 and 8x32 limb layouts where supported.  Capture
hardware counters for cycles, instructions, integer multiply, branches,
cache misses, and vectorization.  Pin cores and report compiler/version/flags.

No CPU result is generalized to GPU/FPGA.

## E13: GPU matrix

For each viable arithmetic/walk combination sweep block size and walkers per
thread.  Record registers/thread, achieved occupancy, local-memory spills,
integer instruction mix, constant/shared/global table traffic, and DP output
traffic.

Run device differential self-tests before performance measurements.  Static
PTX/SASS counts without a device run remain explicitly labelled static.

## E14: FPGA matrix

Compare a generic Montgomery pipeline with a P-256-specific
generalized-Mersenne reducer.  Report:

```
rho steps/s
steps/s/DSP
steps/s/LUT
BRAM usage
Fmax
energy/step (when board power is available)
```

Synthesis estimates are not mixed with measured board throughput.

## E15: factorial ablation and interaction test

Take the individually viable choices and run:

```
baseline
+ arithmetic
+ weak reduction
+ partition
+ table choice
+ batching
+ all combined
```

Then test the most important pairwise interactions, especially:

- fused reduction x weak reduction;
- weak reduction x partition;
- partition x table size;
- batched affine x partition;
- table size x GPU cache/register pressure.

The final table has one row per variant and one cost unit, with baseline and
rho reference included.

## Candidate combined architecture

The high-value combined hypothesis is:

```
many independent affine walkers
+ simultaneous inversion
+ weak P-256 field arithmetic
+ cheap invariant x-derived partition
+ r-adding table
+ sparse distinguished points
```

This is a hypothesis, not a result.  In particular, affine batching is allowed
to lose; the experiment exists to decide it.

## Correctness gates

Every optimized primitive is cross-checked against a slow reference.
Every toy solve must recover the planted logarithm.  Every emitted collision
must verify in the group.  Zero failed relations/collisions are tolerated.
Cycle escapes and restarts are included in cost.

For each final candidate run at least 64 seeded solves per toy size.  Report
confidence intervals because rho solve lengths have high variance.

## Result schema

Every frozen run writes a machine-readable row containing at least:

```json
{
  "curve": "...",
  "n_bits": 0,
  "seed": 0,
  "backend": "...",
  "representation": "...",
  "normalization_K": 0,
  "partition": "...",
  "R": 0,
  "batch_B": 0,
  "dp_bits": 0,
  "fold": 1,
  "verified": false,
  "counted_operations": {},
  "steps": 0,
  "cycle_escapes": 0,
  "restarts": 0,
  "wall_ns": 0,
  "raw_steps_per_s": 0.0,
  "collision_constant": 0.0
}
```

Hardware-specific fields are appended rather than replacing the common
accounting.

## Decision rule

Classify each change using the repository convention:

- **engineering:** same generic boundary, lower complete `S`;
- **relabelling:** a stage gets cheaper while complete `S` does not;
- **accounting:** only the measurement changes;
- **advance:** reserved for a demonstrated movement relative to a mathematical
  floor, not for a faster implementation of rho.

This thread is expected to produce engineering improvements.  It must not call
a constant-factor rho implementation improvement a new ECDLP attack.

## Execution order

Run E1/E2 first because all later experiments consume field arithmetic.
E3/E4 follows, then E5-E8.  E9/E10 is the first architecture-changing
decision.  E11 follows the winning representation.  E12-E14 port the frozen
variants rather than independently changing the algorithm.  E15 is the only
place where combined speedup is claimed.

A negative result is retained with the same seeds and accounting as a positive
one.

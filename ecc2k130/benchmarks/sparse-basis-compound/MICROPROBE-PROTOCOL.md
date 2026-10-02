# Exact-layout sparse-table multicast microprobe: frozen protocol

Source parent: `af613db3bb9e017b69eae830b579748b8e99828c`.

This is a primitive-only admission experiment for the conditional finding in
`RESULTS.md`.  It does not build or time a sparse-basis walk, run a search,
collect distinguished points, solve a collision, or recover a scalar.  The
experiment implementation, field/table construction, correctness checks,
timing, and summarization are native C++17/CUDA.  The existing Modal transport
may copy the clean tree to one GPU and invoke the shell job; it is not part of
the measured interval or scientific computation.

## Frozen question

The earlier static screen transferred 9.2 random-`LDS.U8`
lane-instructions/SM-clock to an unmeasured table pattern.  Measure the exact
pattern instead:

- joint direct-map key `(j,value)`, with `j in 0..7` and a three-bit chunk
  value in `0..7`;
- fixed sparse-to-normal, sparse-to-beta, and beta-to-sparse maps, keyed only
  by the three-bit value; and
- an aligned low-four-word `LDS.128` plus a separate top-word `LDS.32` for
  each of 44 chunks.

The selected field is the statically preferable
`z^131+z^8+z^3+z^2+1`, with beta-root
`0x30a16693fefe59e60962fc4e3ddc388eb`.  All tables are constructed natively at
runtime from that field and the repository's actual beta/normal transforms.
No generated Python artifact is an input.

## Hardware and layout

The job fails closed unless it sees exactly one NVIDIA RTX PRO 6000 Blackwell
Server Edition, compute capability 12.0, 188 SMs, and CUDA compiler 13.3.73
producing native `sm_120` code.

Every timed kernel uses 188 blocks, 512 threads/block, one block/SM, 65,536
rounds, and the route's complete 77,440-byte table allocation plus the
1,024-byte driver reservation.  Each map application executes 44 aligned
128-bit shared loads and 44 top-word shared loads.  Table initialization and
the block barrier are inside the CUDA-event interval; the long loop amortizes
them without silently dropping the real allocation.

The actual arms are:

1. `joint_varying`: actual `L_j`, with lane- and round-varying `(j,value)`;
2. `fixed_normal`: actual sparse-to-normal conversion;
3. `fixed_to_beta`: actual sparse-to-beta conversion;
4. `fixed_from_beta`: actual beta-to-sparse conversion.

Matched controls are:

5. `joint_fixed_j`: the same joint table and instructions with `j=0` for all
   lanes while values still vary;
6. `joint_uniform_key`: the same joint table and instructions with one
   `(j,value)` broadcast by every lane in a warp; and
7. `fixed_normal_clone`: the exact `fixed_normal` kernel and inputs, providing
   the A/A noise pair.

The controls diagnose joint-key bank conflicts and multicast benefit.  They
do not enter the sparse-route capacity calculation.

## Correctness and native binding

Timing is suppressed unless all of the following pass:

1. Native field construction rechecks both conversion ranks, every map rank,
   all 131 basis directions, and 4,096 deterministic dense inputs against the
   repository arithmetic.
2. GPU one-application outputs match the independent host result for every
   basis and dense input for the joint and three fixed maps.
3. A GPU key histogram for the varying joint arm observes all 64 `(j,value)`
   keys outside the timed interval.
4. Every timed arm stores all thread outputs.  One spread-selected thread per
   arm is replayed natively through all 65,536 rounds, and every later sample
   must reproduce the first sample's complete-output SHA-256.
5. Runtime attributes report zero local bytes, at least one resident
   512-thread block/SM, and exactly 77,440 dynamic shared bytes requested.
6. The retained native binary, build command, ptxas log, cubin/SASS, source,
   tables, samples, and result are SHA-256 bound.  A native SASS auditor must
   find both `LDS.128` and scalar `LDS` in every table kernel, reject local
   loads/stores and calls in those kernels, and bind its report to the exact
   executable/cubin hashes.

The benchmark evolves each input between maps and feeds the table output into
the next round so the loads cannot be hoisted or deleted.  The same evolution
and output folding are used in every arm.

## Timing and decision

One warmup per arm is excluded.  Five measured rounds alternate forward and
reverse arm order.  `fixed_normal` and `fixed_normal_clone` swap pair order on
successive rounds.  Each row must report the exact applications, 88 LDS
instructions/application, finite positive CUDA-event time, zero errors, and
the bound output digest.

Let `r_j`, `r_n`, `r_b`, and `r_s` be the medians in map applications/s for
`joint_varying`, `fixed_normal`, `fixed_to_beta`, and `fixed_from_beta`.
The exact lookup-only sparse-route capacity is

```
R_lookup = 1 / (2/r_j + 1/r_n + (1/16)*(1/r_b + 1/r_s)).
```

The primitive passes the prerequisite only when:

- all correctness, binary, SASS, resource, and digest gates pass;
- maximum absolute A/A pair drift is below 1%; and
- `R_lookup >= 1.05 * 15.436677 B/s = 16.20851085 B/s`.

A pass justifies a separately frozen whole-walk prototype.  A failure retains
the current conditional no-go.  Neither outcome is a whole-walk throughput
measurement.

## Frozen `PACKED_ALU_SQUARE=1` price

The current compound ledger has 38.125 `CLMAD` lane-instructions/update:
510 from 85 field products, 20 from five normal-basis inverse squarings, and
80 from sixteen five-word polynomial lambda squarings, all per batch 16.

Moving each lambda square to the sparse ALU spread removes exactly five
`CLMAD`s/update:

```
38.125 - 5 = 33.125 CLMAD/update.
```

At 188 SMs, 2.430 GHz, and the ideal two-`CLMAD` lane rate, its isolated
ceiling is **27.582792 B/s**, 1.0609 times the 26 B/s objective.  The native
source circuit adds five `spread32alu` calls, each five
shift/OR/AND stages, or **75 source ALU operations/update**, before compiler
fusion.  The previously charged sparse reducer and maps therefore move from
`476.625 + 990.625 = 1,467.25` to **1,542.25 source operations/update** in
those stages, still below the reference subledger's 2,276.25.  These are
source operations, not SASS or a rate prediction.

This arithmetic price is frozen before the lookup measurement.  A lookup
pass does not prove 26 B/s: a later whole-walk candidate must compile without
spills and jointly fit the measured lookup, ALU, `CLMAD`, state, and control
costs.

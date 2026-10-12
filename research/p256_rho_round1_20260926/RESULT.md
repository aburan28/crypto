# P-256 rho arithmetic round 1 — pilot result

**Status:** bounded pilot.  This is not a P-256 discrete-log solve and not an
end-to-end P-256 speedup claim.

The pre-registered program is #836 / #852.  This round adds an executable
harness and runs the parts that can be measured without a CUDA device or FPGA
toolchain.

## What was run

Full-width P-256 field multiplication was checked on 10,000 deterministic
random pairs using

```
p = 2^256 - 2^224 + 2^192 + 2^96 - 1
seed = 0x50323536
```

The generalized-Mersenne fold was compared exactly with Python big-integer
`(a*b) % p` on every pair.

For collision-quality testing, complete Pollard-rho solves were run on the
prime-order toy curve

```
p = 8191
E: y^2 = x^3 + x + 5
#E = 8053 (prime)
```

with 16 deterministic planted secrets per cell.  The walk used affine points,
the negation map, distinguished points, and the same r-adding table construction
within each matched cell.

## Field result (stage diagnostic)

| backend | median ns/mul on this Python host | correctness |
|:--|--:|:--|
| Python big-int `a*b % p` | 484.8 | reference |
| explicit Solinas fold | 2780.8 | 10,000/10,000 |

The explicit Python fold is **5.74x slower** than Python's native modular
reduction.  This says nothing about the CUDA/C++ hypothesis: Python big
integers execute optimized native code while the explicit fold executes a
Python loop.  It is retained as a negative stage diagnostic and correctness
receipt, not used to select E1/E2.

## Complete toy rho panel

Every one of the 192 planted DLPs verified.

| R | partition | verified | mean counted walk steps |
|--:|:--|--:|--:|
| 8 | xlow | 16/16 | 989.1 |
| 8 | xor2 | 16/16 | 1657.2 |
| 8 | mix | 16/16 | 1010.1 |
| 16 | xlow | 16/16 | 724.6 |
| 16 | xor2 | 16/16 | 572.1 |
| 16 | mix | 16/16 | 581.2 |
| 32 | xlow | 16/16 | 414.0 |
| 32 | xor2 | 16/16 | 376.0 |
| 32 | mix | 16/16 | 521.9 |
| 64 | xlow | 16/16 | 263.9 |
| 64 | xor2 | 16/16 | **263.1** |
| 64 | mix | 16/16 | 287.4 |

This pilot therefore says only that, in this small implementation, increasing
the table from 8 to 64 removes a large amount of walk inefficiency and that
`xlow` and `xor2` are indistinguishable at the best tested table size within
the noise of 16 rho solves.  It does **not** establish that xor-mixing is better
than the essentially-free low-limb partition.

## E1–E15 status after this round

| experiment | round-1 status |
|:--|:--|
| E1/E2 field backends + fused Solinas | correctness oracle + Python diagnostic; native CUDA/C++ comparison still open |
| E3/E4 weak reduction | harness implemented; native overflow/instruction accounting still open |
| E5 partition cost | harness implemented; native instruction measurement open |
| E6 partition quality | pilot executed |
| E7 representation invariance | affine pilot is canonical; Jacobian rescaling gate remains required |
| E8 table sweep | pilot executed, R=8..64 |
| E9/E10 Jacobian/batched affine | repository already has simultaneous-inversion rho machinery; P-256-specific matched run remains open |
| E11 DP batching | parameterized in harness; full density matrix remains open |
| E12 CPU | Python host pilot only; native 4x64/8x32 matrix open |
| E13 GPU | **not measured** — requires CUDA device |
| E14 FPGA | **not measured** — requires synthesis/board toolchain |
| E15 factorial ablation | pilot table/partition interaction only; final ablation remains open |

## Decision

No optimization is promoted from this pilot.  The next decisive measurement is
native E1/E2 in `gpu/ecc`: specialized Montgomery versus eager/fused P-256
Solinas under the actual rho point-add sequence.  E8's small-group result makes
`R=32,64,128,256` the useful native table range; R=8 is retained as a control.

**Classification:** engineering diagnostic.  The generic-group boundary did
not move.

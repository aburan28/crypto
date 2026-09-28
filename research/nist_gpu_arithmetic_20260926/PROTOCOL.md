# Protocol: RTX PRO 6000 P-256/P-384 arithmetic tournament

Status: pre-registered before any RTX PRO 6000 measurement from this branch. Amendment 1 (still pre-measurement) adds the hand-written sparse Montgomery reduction carry chain to candidate 2.

## Question

Which arithmetic representation and launch configuration maximizes useful
P-256 and P-384 field/point throughput on an NVIDIA RTX PRO 6000 Blackwell
(sm_120), without trading correctness for an isolated microbenchmark?

## Frozen hardware/software target

- GPU: RTX PRO 6000 Blackwell, compute capability 12.0.
- CUDA image: 13.3.1 devel Ubuntu 24.04.
- 32-bit limbs only in this round.
- deterministic input seed: 0x50323536.
- field batch: 2^20 elements.
- P-256: 128 field iterations per timed launch.
- P-384: 64 field iterations per timed launch.
- point iterations: field iteration count / 16.

NVIDIA documents 64K 32-bit registers per sm_120 SM and at most 48 resident
warps. Register pressure is therefore part of the result, not a post-hoc
excuse.

## Frozen candidates

Arithmetic:
1. sparse CIOS Montgomery, portable product row;
2. sparse CIOS Montgomery, explicit PTX product and sparse-reduction rows;
3. direct generalized-Mersenne/Solinas, portable product row;
4. direct generalized-Mersenne/Solinas, explicit PTX product row.

The sparse CIOS reduction uses n'=1 and rewrites every modulus-row multiply
using the signed radix-2^32 limb in {0,+1,-1,-2}. The direct Solinas candidate
uses the published prime identities but the same product-row primitive.

GPU configurations:
- threads/block = 64, 96, 128, 160, 192, 256;
- portable: unrestricted registers;
- PTX: max registers = unrestricted, 96, 128, 160;
- three repeats per cell, median selected.

Workloads:
- dependent field multiply (latency);
- two independent multiply chains (ILP/throughput);
- dependent squaring;
- a=-3 Jacobian doubling;
- mixed Jacobian/affine addition.

## Correctness gates

Before device timing:
- 20,000 deterministic random field pairs per curve and representation must
  match Boost.Multiprecision for mul, square, add, subtract, and conversion;
- P-256 and P-384 Jacobian doubling and mixed-add formulas must match an
  independent affine big-integer calculation.

The device build must run successfully for both portable and PTX paths before
a PTX result can be selected.

Any correctness mismatch rejects the configuration regardless of throughput.

## Selection rule

For each curve, keep separate winners for:
- field multiply throughput;
- point doubling;
- mixed addition.

The point-arithmetic winner is primary. An isolated field winner is not called
the curve winner if register pressure makes point arithmetic slower.

A speed claim requires the median of three repeats and at least a 5% advantage
over the next simpler accepted configuration. Smaller differences are ties and
prefer the simpler implementation.

## Follow-up, not part of this outcome

After the first winner is frozen:
- Nsight Compute counters for integer issue, eligible warps, occupancy and
  local-memory spills;
- 4-chain ILP if register headroom remains;
- dedicated squaring only if SASS shows duplicated products survive;
- fuse field operations across point formulas;
- compare batch-affine inversion only in a complete ECC workload.

No result from this protocol is an ECDLP complexity claim.

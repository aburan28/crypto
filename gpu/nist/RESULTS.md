# gpu/nist first device run — RTX PRO 6000 Blackwell

One NVIDIA RTX PRO 6000 Blackwell Server Edition (cc 12.0, 188 SMs),
CUDA 13.3, `-arch=sm_120`, through `gpu/nist/modal_app.py` (Modal app
`nist-prime-gpu-arithmetic`; runs linked below). Two device builds from the
same source: `NIST_PTX=1` (mad.lo/hi.cc carry chains in `nist_ptx.cuh`)
and `NIST_PTX=0` (portable software product, the control). CPU suites
(`make test`) pass in the same image before every device run.

**What this round fixed to get a first honest number.** The kernel had never
run on a GPU, and two compile-time/runtime bugs made the first run measure
nothing:

1. `F::batch_inv` was host-only (`static`, no `NF_HD`). nvcc warns
   `calling a __host__ function ... is not allowed` and drops the call in
   device code, so `rho_batch_step` compiled to an empty body: `k_rho`
   was 4 registers / 0 stack and ran in microseconds with 0 DPs. The CPU
   test never saw it because g++ compiles host-only. `cuobjdump -sass`
   showed `LDC R1; EXIT`. Fix: `NF_HD` on `batch_inv`
   (`nist_fields.cuh:132`).
2. `rho_count`'s device branch had no null check; the bench passes
   `cycles = nullptr` and fruitless 2-cycles occur at ~1/(fold·R) ≈ 1/512
   per step, so the first cycle escape faulted with an illegal memory
   access. Fix: `if(ctr)` guard, matching the host branch.

**Correctness evidence on hardware.** Every rho cell reports its DP count;
the walk is a deterministic function of the point and DPs occur with
probability 2^-24 per step. Over the sweep, observed/expected DPs were
e.g. 11,969/12,085 (P-256 first run) and 6,094/6,016, 5,966/6,016 (best
cells) — all within binomial noise. The CPU suite additionally cross-checks
the new compile-time-W stepper `rho_batch_step_w<Wc>` bit-for-bit against
the generic stepper for W = 2, 4, 8, 16 (`test_rho`, 4/4 identical).

## Field arithmetic (threads=128, first run)

| | P-256 Montgomery | P-256 Solinas | P-384 Montgomery | P-384 Solinas |
|---|---:|---:|---:|---:|
| mul chain1, Gmul/s (PTX=1) | 73.4 | 70.4 | 32.2 | 36.7 |
| mul chain1, Gmul/s (control) | 54.9 | 61.2 | 20.3 | 30.7 |
| point_double, M/s (PTX=1) | 6084 | 6550 | 1712 | 1647 |

## Rho walk

First working measurement (runtime-W kernel, threads=128, W=8):

| curve | control Gstep/s | PTX Gstep/s | PTX ratio |
|---|---:|---:|---:|
| P-256 | 1.075 | 1.282 | 1.19 |
| P-384 | 0.210 | 0.289 | 1.38 |

Compile-time-W stepper (`rho_batch_step_w<Wc>`, `k_rho<C,Wc>`), PTX build,
threads/block × W sweep (8192 iters per cell):

| threads | W | P-256 Gstep/s | P-384 Gstep/s |
|---:|---:|---:|---:|
| 64 | 2 | 0.454 | 0.091 |
| 64 | 4 | 0.772 | 0.162 |
| 64 | 8 | 1.194 | 0.265 |
| 64 | 16 | 1.271 | 0.410 |
| 128 | 2 | 0.464 | 0.095 |
| 128 | 4 | 0.786 | 0.169 |
| 128 | 8 | 1.217 | 0.279 |
| 128 | 16 | 1.306 | 0.434 |
| 256 | 2 | 0.466 | 0.096 |
| 256 | 4 | 0.812 | 0.171 |
| 256 | 8 | 1.24 | 0.32 |
| 256 | 16 | **1.340** | **0.438** |

k_rho instantiations (PTX build): P-256 W=2/4/8 at 92/116/116 registers,
W=16 at 180; P-384 198–204 registers; stack 688–2576 bytes growing with W.

**Verdict (engineering, per the boundary rule): the walk got faster without
changing the algorithm.** From the first working measurement (control
build, threads=128, W=8) to the tuned cell (PTX, threads=256, W=16):
P-256 1.075 → 1.340 Gstep/s (**1.25×**) and P-384 0.210 → 0.438 Gstep/s
(**2.09×**). Step count is unchanged; S = steps/√n is flat.

## Round 2: register-resident inversion, W=32, maxrregcount

`rho_batch_step_w` now passes the den/scratch arrays of the batched
inversion by reference (`F::batch_inv_w<Wc>`, unrolled; `test_rho` checks
it against `batch_inv` on 2000 random inputs per W), and `k_rho<C,Wc>`
gained a W=32 instantiation. Swept threads × W unrestricted, then
maxrregcount 96/128/160 (PTX build, 8192 iters per cell):

| threads | W | P-256 Gstep/s | P-384 Gstep/s |
|---:|---:|---:|---:|
| 128 | 8 | 1.212 | 0.279 |
| 128 | 16 | 1.310 | 0.408 |
| 128 | 32 | 1.603 | 0.538 |
| 256 | 16 | 1.339 | 0.417 |
| 256 | 32 | **1.632** | **0.548** |
| 256, maxreg=96 | 32 | 1.589 | 0.447 |
| 256, maxreg=96 | 16 | 1.601 | 0.448 |
| 128, maxreg=128 | 32 | 1.438 | 0.418 |
| 256, maxreg=128 | 32 | 1.595 | 0.455 |
| 256, maxreg=160 | 32 | 1.595 | 0.455 |

Every register cap loses to the unrestricted build (96/128/160 trade spill
for occupancy and come out behind), so the 255-register W=32 instantiation
is the keeper. DPs match the 2^-24 rate in every cell (best cell: 12,154
observed / 12,096 expected).

**Attribution is honest but not separated**: round 1 never ran W=32, so the
round-2 gain cannot be cleanly split between the larger batch and the
reference-passed inversion. What is separated: at fixed W ≤ 16 the two
rounds agree within noise (128thr W=16: 1.306/0.434 round 1 vs
1.310/0.408 round 2), and every maxreg cap regresses. The W=32 cell wins
on both curves.

Cumulative from the first working measurement (control, threads=128,
W=8) to the current best (PTX, threads=256, W=32, unrestricted):
P-256 1.075 → 1.632 Gstep/s (**1.52×**) and P-384 0.210 → 0.548 Gstep/s
(**2.61×**). Still engineering: the walk computes the same function on the
same curves, one step per counted step, S flat.

Next rungs: W=64 (the 255-register ceiling now binds: W=32 sits at it for
both curves, so W=64 will spill and probably regress — worth one cell to
confirm), in-walk reseeding off the critical path, and a second walk per
step from the same table read.

End-to-end scale (practicality note, not a claim): at 1.63 Gstep/s a
P-256 collision needs ~√(π·2^256/4)/1.63e9 s ≈ 1.9e30 s, so these rates are
attack-layer throughput measurements on retired-practice targets, not a
statement about NIST curves being weak.

Runs (Modal): first device run ap-cn0g0oMWh2SBtpWS5RrCEf (empty-kernel
discovery), control build (validate), PTX build ap-39mdjitYFxkMRHk6pFdd2P,
round-1 sweep ap-o766pq9fqlurbzSYcqZoe6, round-2 sweep
ap-0KNqNM8D1J4YQn5z67Ucc0 (timed out at Modal's 4 h during the last
maxreg=160 cells; all cells reported above completed normally).

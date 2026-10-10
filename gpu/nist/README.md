# gpu/nist — P-256 / P-384 arithmetic on RTX PRO 6000 Blackwell

This directory is a clean arithmetic target for NIST P-256 and P-384. It is separate from gpu/ecc because that code is fixed at 256-bit limbs and grew around secp256k1's prime and rho walk.

## Architecture

RTX PRO 6000 Blackwell is an sm_120 device. The design uses 32-bit limbs (8 for P-256, 12 for P-384), thread-private field elements, CIOS Montgomery multiplication, and a=-3 Jacobian formulas.

Both primes have low limb 0xffffffff, hence Montgomery n' = 1. More importantly, in radix B=2^32 every modulus limb is 0, 1, B-1, or B-2. The reduction row therefore uses carry arithmetic instead of integer multiplication. P-256's signed limb pattern is [-1,-1,-1,0,0,0,+1,-1]; P-384's is [-1,0,0,-1,-2,-1,-1,-1,-1,-1,-1,-1].

The CUDA path adds explicit mad.lo.cc / madc.hi.cc product-row carry chains and hand-written sparse Montgomery reduction carry chains. A direct generalized-Mersenne/Solinas field is retained as a matched competitor, so the RTX PRO measurement decides rather than assuming either representation wins.

## Pollard-rho ECDLP walk

On top of the field layer this directory now carries a full Pollard-rho engine for the NIST prime curves, so the arithmetic target is also an ECDLP solver that matches the CPU engine in `src/cryptanalysis/ecdlp_nist/`.

| File | What it is |
|---|---|
| `nist_curve.cuh` | Affine points and the small curve arithmetic the walk needs, generic over the field with a pluggable coefficient `a` (`AMinus3` for P-256/P-384). No Jacobian and no a=-3 assumption in the walk layer: the hot step never inverts. |
| `nist_rho.cuh` | The r-adding walk: negation-map fold (fold 2, `sqrt(2)` fewer steps; j != 0 so `+-1` is all of Aut), distinguished points, fruitless-cycle detection and escape, splitmix64 seeds, SoA batched state, and `rho_batch_step` batching the one inversion across a thread's walks. Generic over `F::N`, so one body serves both curves. |
| `nist_rho_host.hpp` | Host driver: jump-table build, the walk carrying each walk's `(a,b)` with `P=a*G+b*Q` in lockstep with `phase_b`, and collision -> logarithm. |
| `test_rho.cpp` | CPU verification (below). |

The CUDA kernels (`k_rho_seed`, `k_rho`) and their runner live in `bench.cu` and call the same `rho_batch_step` the CPU test exercises.

## Correctness

Run: cd gpu/nist && make test

`test_cpu` checks 20,000 random field pairs on each curve against Boost.Multiprecision and checks Jacobian doubling and mixed addition against independent affine arithmetic (Boost is a CI dependency).

`test_rho` verifies the curve/walk/solver headers with **no external dependency**, on group-law identities that need no big-integer library: `a*inv(a)=1` and `batch_inv==inv`; `P+O=P`, `P+(-P)=O`, `[2]P` via add, associativity, `[a]G+[b]G=[a+b]G`, and `affine_add_with_inv==add`; the negation-map fold is a class function with the right aut code; the single step is deterministic; `rho_batch_step` (the exact kernel body) runs on the toy and emits only valid on-curve distinguished points; and the decisive end-to-end DLP solve on a ~2^20 prime-order toy curve whose `(n,G,Q,k)` were checked independently at construction. The folded walk solves it in fewer steps than the unfolded one, the expected `sqrt(2)`.

**No kernel here has been run on a GPU** (the environment has no NVIDIA device). `make test` compiles the same `.cuh` headers the kernel uses with g++, so a green run means the field/curve/walk/solver logic is correct; only the launch configuration and memory layout remain GPU-specific. `make bench` needs `nvcc`.

## RTX PRO 6000 benchmark

Run: make bench ARCH=sm_120 NIST_PTX=1 && ./bench 128

The benchmark reports dependent and two-chain field multiplication, squaring, a=-3 doubling and mixed addition for both Montgomery and Solinas fields. The tuner sweeps 64..256 threads/block and PTX register caps 96/128/160 plus unrestricted. sm_120 has 64K 32-bit registers per SM and at most 48 resident warps, so especially for P-384 the winner cannot be selected from instruction count alone.

The Modal runner automates portable-vs-PTX and block-size tuning: modal run modal_app.py::validate, then modal run modal_app.py::tune. NIST_GPU defaults to RTX-PRO-6000.

## Measurement rule

Do not call the highest isolated Gmul/s configuration the point-arithmetic winner. Keep separate rows for field multiply, point double and mixed add, and carry ptxas register counts beside throughput. The next round should add Nsight Compute counters for integer issue, eligible warps, local spills and occupancy, then test 1/2/4 independent arithmetic chains per thread.

## Sources

NVIDIA RTX PRO 6000 Blackwell Server Edition: https://www.nvidia.com/en-us/data-center/rtx-pro-6000-blackwell-server-edition/

NVIDIA Blackwell Tuning Guide: https://docs.nvidia.com/cuda/blackwell-tuning-guide/

CUDA Programming Guide: https://docs.nvidia.com/cuda/cuda-programming-guide/

# gpu/nist — P-256 / P-384 arithmetic on RTX PRO 6000 Blackwell

This directory is a clean arithmetic target for NIST P-256 and P-384. It is separate from gpu/ecc because that code is fixed at 256-bit limbs and grew around secp256k1's prime and rho walk.

## Architecture

RTX PRO 6000 Blackwell is an sm_120 device. The design uses 32-bit limbs (8 for P-256, 12 for P-384), thread-private field elements, CIOS Montgomery multiplication, and a=-3 Jacobian formulas.

Both primes have low limb 0xffffffff, hence Montgomery n' = 1. More importantly, in radix B=2^32 every modulus limb is 0, 1, B-1, or B-2. The reduction row therefore uses carry arithmetic instead of integer multiplication. P-256's signed limb pattern is [-1,-1,-1,0,0,0,+1,-1]; P-384's is [-1,0,0,-1,-2,-1,-1,-1,-1,-1,-1,-1].

The CUDA path adds explicit mad.lo.cc / madc.hi.cc product-row carry chains and hand-written sparse Montgomery reduction carry chains. A direct generalized-Mersenne/Solinas field is retained as a matched competitor, so the RTX PRO measurement decides rather than assuming either representation wins.

## Correctness

Run: cd gpu/nist && make test

The test checks 20,000 random field pairs on each curve against Boost.Multiprecision and checks Jacobian doubling and mixed addition against independent affine arithmetic. The host path always uses portable uint64_t, so it remains the reference for the device specialization.

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

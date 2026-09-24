# gpu/fes — Gray-code FES enumeration on the GPU

Fast exhaustive search (FES, Bouillaguet et al.; libfes-lite) enumerates every
assignment of a quadratic Boolean system in Gray-code order, so consecutive
points differ in one variable and each equation is updated from a small table
of derivatives. In a binary index-calculus decomposition (Weil descent /
summation polynomials over `F_2`), FES is the **relation-search inner loop**,
and it is embarrassingly parallel: the `2^n` cube of assignments splits into
independent sub-cubes, one per GPU thread.

This directory is that kernel, in three forms that share one algorithm:

| file | what it is |
|---|---|
| `fes.cuh` | the algorithm: packed quadratic evaluation, the derivative-maintained Gray walk, and `fes_search_subcube` — one sub-cube (a fixed assignment of the top `n-low_bits` variables). Compiles as host C++ **and** as CUDA. |
| `fes.cu` | CUDA launcher: one thread per sub-cube, results appended through an atomic counter. |
| `fes.metal` | Metal compute shader: the same kernel for Apple silicon. |
| `test_cpu.cpp` | host verification — runs `fes_search_subcube` for every prefix and checks the collected solutions against a from-scratch brute force. **No GPU required.** |
| `fref.py` | reference generator/oracle: emits a random system with a planted solution and its independently brute-forced solution set. |

## Build & test

```sh
make test                 # host verification, no GPU (the correctness gate)
make test N=20 M=22 SEED=3
make cuda                 # compile+link the CUDA launcher (needs nvcc)
make nvcc-check           # compile-only for sm_70/80/90
make metal                # compile the Metal shader (macOS + xcrun)
```

`make test` compiles the exact `fes.cuh` the device kernels use and proves the
field-free evaluation, the Gray-code derivative walk, and the sub-cube sharding
(it checks that 1, 4, 16 and 64 shards all find the same solution set as brute
force). Only the launch is device-specific.

## Honest scope

Like the rest of `gpu/` in this repository, these kernels are **compiled and
cross-checked, not yet wired into the Rust `icx`/`ic` binaries**. The Rust
pipeline's FES solver (`crypto_lib::cryptanalysis::mq_fes`, the `fes-f2` /
`fes-f2-wide` plug-ins) is the production path today; it already uses AVX2 /
AVX-512 on x86. This directory is the offload target: the CPU FES and this
kernel agree by construction (same algorithm), and wiring is a file-exchange
step (emit the packed system, launch, read back candidate solutions and
re-verify them in the group) matching how `gpu/ecc2k` exchanges its pair table.

No GPU was available in the environment that wrote this; correctness rests on
the host verification and the CUDA/Metal **compile** checks in CI. A real-GPU
run is future work, as it is for the other kernels here.

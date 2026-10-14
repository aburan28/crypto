# gpu/rootfind — CUDA root-finding 2-decomposition oracle (Blackwell sm_120)

Batched `R = ε·F_i + τ·F_j` decomposition over a factor base by the
Semaev-`S₃` root identity (`x(R ∓ F_i)` lookup), for a 64-bit prime
field mirroring `prime_fast::Fp64` and a 4×64-limb P-256 field mirroring
`ecc::p256_field`. Design, validation results, throughput **estimates**
and the claim boundary: [`RESEARCH_GPU_ROOTFIND_DECOMPOSITION.md`](../../RESEARCH_GPU_ROOTFIND_DECOMPOSITION.md).

**Not run on a GPU.** Written on a machine without an NVIDIA device or
`nvcc`. The kernels' per-thread function is `__host__ __device__` and is
what the host emulation executes; that emulation is validated bit-for-bit
against the CPU `prime_fast` oracle (`rs/src/bin/gpu_rootfind_validate.rs`).

| File | What |
|---|---|
| `fp64.cuh` | Montgomery `F_p`, `p < 2^62`, one `u64` (= `Fp64`) |
| `fp256x64.cuh` | P-256 Montgomery on 4×`u64`, SOS (= `P256FieldElement`) |
| `rootfind.cuh` | per-thread oracle, `SortedTable` / `HashTable`, kernel args |
| `rootfind_host.hpp` | table construction, conversions, emulation driver |
| `rootfind_api.h/.cpp` | C ABI (`rf64_decompose`, `rf256_decompose`, …) |
| `rootfind_cuda.cu` | CUDA kernels + device driver (sm_120 notes in header) |
| `bench.cu` | `selftest` / `throughput` / `p256` on a device |
| `test_cpu.cpp` | host self-test (`make test`) |
| `Makefile`, `CMakeLists.txt` | CUDA optional in both |
| `rs/` | separate cargo crate `gpu-rootfind`: raw FFI, safe wrappers, `build.rs`, `gpu_rootfind_validate` harness (depends on the root crate; the root crate has no feature/dep/build script for it) |

```bash
make test                      # no GPU needed
make bench ARCH=sm_120 && ./bench selftest     # on the RTX PRO 6000 (untested here)
cargo test --release --manifest-path rs/Cargo.toml          # FFI + emulation vs CPU reference
cargo run --release --manifest-path rs/Cargo.toml --bin gpu_rootfind_validate
```

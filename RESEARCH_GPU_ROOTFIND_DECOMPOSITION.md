# Research Note: GPU (CUDA) Root-Finding 2-Decomposition Oracle for Blackwell sm_120

**Date.** 2026-10-12. **Status.** Host-emulated kernel validated fail-closed
against the CPU `prime_fast` oracle on public synthetic known-answer
instances (§4). **The CUDA path has NOT been executed on a GPU**: this was
written on an Apple-silicon Mac with no NVIDIA device and no `nvcc`
(`which nvcc` → not found). Every throughput number for sm_120 in §3 is an
**ESTIMATE** from published hardware parameters, not a measurement.
**No security claim** (§6).
**Code.** `gpu/rootfind/` (headers shared by device and host, CUDA backend,
C ABI, CPU self-test, Makefile + CMake), `gpu/rootfind/rs/` (separate
cargo crate `gpu-rootfind`: FFI, safe wrappers, build script, and the
`gpu_rootfind_validate` harness), `src/cryptanalysis/gpu_rootfind.rs`
(root crate: `Relation` type, CPU reference, group re-add check);
evidence under `research/gpu_rootfind_decomposition_20261012/`.
**Builds on.** `RESEARCH_SEMAEV_ROOTFIND_VS_GROEBNER.md` (the identity
`roots_X S₃(x_R, x_i, X) = {x(R+F_i), x(R−F_i)}`, Lean-checked) and
`prime_fast::collect_relations`.

---

## 0. Summary

The root-finding decomposition oracle is embarrassingly parallel: for a
target `R` and factor base `F = {F_0 … F_{N−1}}` it is `N` independent
`(target, i)` cells, each two slope computations, one squaring and two
table probes, with the only shared work being one field inversion per
target (amortised by Montgomery's trick). This note ports it to CUDA for
the NVIDIA RTX PRO 6000 Blackwell (sm_120) and, because no GPU was
available, validates it through a **host emulation that executes the
kernel's per-thread function unchanged** (`__host__ __device__` headers,
the convention of `gpu/ecc/`).

- **64-bit field** (`p < 2^62`, mirroring `prime_fast::Fp64` constant for
  constant) and a **4×64-limb P-256 field** (mirroring
  `ecc::p256_field` / `ct_bignum::mont_mul`). Both kernels are
  instantiations of one templated per-thread function.
- **Validation (emulation, this machine):** on `a = −3` prime-order curves
  at 40, 56 and 61 bits with factor bases of 4 096–16 384 points, over
  **21 811 targets in the three evidence files** (1 280 positive controls
  `±F_i ± F_j`, 48 degenerate controls `±F_i`, 3 production targets from
  `collect_relations`, 20 480 random on-curve rejections; ≈ 24 600 with
  the 48-bit console run and the unit tests), all **6 (lookup × chunk)
  configurations** returned exactly the CPU `prime_fast` relation set —
  **1 331 relations, 0 missing, 0 spurious, 0 duplicates, 0 re-addition
  failures, 0 controls missing** — and every relation
  re-adds to its target in the group. The P-256 instantiation reproduces
  positive controls built with the crate's `P256FieldElement`.
- **Estimate for sm_120 (§3):** ≈ 5×10¹⁰–10¹¹ `(target, i)` cells/s for
  the 64-bit field, i.e. ~10⁶ targets/s at `|F| = 2^16`; roughly
  300–1 000× the multi-threaded host emulation measured here
  (1.0–1.9×10⁸ cells/s). To be falsified by `./bench throughput` on the
  device.

## 1. Design

### 1.1 Mathematics (unchanged)

For `y² = x³ + ax + b` and a target `R = (x_R, y_R)`, `R = ε F_i + τ F_j`
with `ε, τ ∈ {±1}` iff `x(R − ε F_i) = x_j`. With `d = x_R − x_i`:

```
λ₊ = (y_R + y_i)·d⁻¹      ε = +1   (computes R − F_i)
λ₋ = (y_R − y_i)·d⁻¹      ε = −1   (computes R + F_i)
x₃ = λ² − x_R − x_i ;  hit ⇔ x₃ ∈ {x_j}
y₃ = λ(x_R − x₃) − y_R ;  τ = +1 if y₃ = y_j else −1
```

This is the loop body of `prime_fast::collect_relations`; the kernel
emits `(k, i, ε, j, τ)` for `j ≥ i` and `(k, i, ε, i, 0)` for the
degenerate `R = ε F_i`, the same convention.

### 1.2 Thread mapping and memory layout (`rootfind.cuh`)

- **One thread = one target × one `CHUNK`-long slice of F.** Thread
  `t → (k = t mod M, slice = t div M)`. Consecutive lanes hold
  consecutive *targets* and the *same* slice: the target loads are
  coalesced (128 threads read 1 KiB of contiguous `(x_R, y_R)` for Fp64)
  and every factor-base load inside the slice loop is **warp-uniform** —
  one broadcast per warp rather than 32 scattered loads.
- **SoA factor base**: `x[]` and `y[]` separate; a probe miss (the common
  case) never touches `y`.
- **Inversion: Montgomery's trick per thread over its slice.** One Fermat
  inversion `a^(p−2)` (fixed public exponent → no warp divergence; ~62
  squarings + ~31 multiplications at 62 bits) plus `3(CHUNK−1)`
  multiplications, against ~95 multiplications per cell for one Fermat
  inversion per `(target, i)`. At `CHUNK = 32` that is ≈ 6 vs 95
  multiplications per cell. The trade-off is register pressure: the
  thread holds `2·CHUNK` field elements (`d[]` and the prefix products),
  i.e. `4·CHUNK` 32-bit registers for Fp64 and `16·CHUNK` for P-256,
  which is what sets occupancy (§3). Binary extended Euclid, which the
  Rust side uses, was rejected for the device: data-dependent trip counts
  serialise a warp. Both give the same field element, so the Montgomery
  representation (unique) matches `Fp64::inv` bit for bit.
- **Lookup keyed by the Montgomery form itself** (its low limb for
  P-256), never the canonical value: the representative is unique, so
  equality of representations is equality of values, and the
  `from_mont` multiplication per probe is saved. A 64-bit key collision
  (possible only for P-256) is rejected by a full-width compare.
- **Two table shapes behind one probe protocol:**
  `SortedTable` (ascending keys + index, lower-bound binary search,
  `log₂N` dependent loads, global memory, any `N`) and `HashTable`
  (open addressing, linear probing, load ≤ ½, ~1–2 loads). The CUDA
  kernel stages the hash table in **static shared memory** when it fits
  `RF_SHARED_HASH_CAP = 4096` slots (48 KiB ⇒ `|F| ≤ 2048`); above that
  the default (`RF_LOOKUP_AUTO`) is the global sorted table. The
  emulation applies the same rule so both backends run the same table.
- **Output**: relations appended through one global `atomicAdd`. Hits
  are rare (≈ `2N²/p` per random target) so it is not a contention
  point. If the buffer fills the counter keeps counting; the host sees
  `n_found > cap` and the call fails closed (`RF_ERR_OVERFLOW`); the Rust
  wrapper retries once with the exact count.

### 1.3 Fields

| Header | Mirrors | Representation | Device/host split |
|---|---|---|---|
| `fp64.cuh` | `prime_fast::Fp64` (`R = 2^64`, `p < 2^62`) | one `u64`, Montgomery | `__umul64hi` vs `unsigned __int128` |
| `fp256x64.cuh` | `ecc::p256_field` + `ct_bignum::mont_mul` (SOS, `R = 2^256`) | 4×`u64`, Montgomery | same |

Constants (`R mod p`, `R² mod p`, `−p⁻¹ mod 2^64`) are either computed
exactly as the Rust code does (`fp64_ctx_new` ≡ `Fp64::new`) or copied
from `p256_field.rs` and re-derived in `test_cpu.cpp` from the
definition. The P-256 header uses 64-bit limbs **to mirror the crate
field, not for speed**: on NVIDIA hardware 32-bit CIOS with inline-PTX
carry chains (`gpu/ecc/fp256.cuh`, 352 PTX instructions per multiply)
is the production choice, and the kernel is templated on the field type
so that it can be dropped in (§5, future work).

### 1.4 Host glue and FFI

- `rootfind_api.h/.cpp` — C ABI: `rf64_decompose`, `rf256_decompose`,
  `rf_cuda_available`, `rf_backend_info`, `rf_strerror`. Canonical
  coordinates in, index relations out. Backend `AUTO|EMULATE|CUDA`,
  lookup `AUTO|SORTED|HASH`, `chunk ∈ {16, 32, 64}`. Rejects even/out-of-
  range `p`, non-canonical coordinates, bad chunk; reports overflow.
- `rootfind_cuda.cu` — the three kernels (`k_decompose_sorted`,
  `k_decompose_shared_hash`, `k_decompose_global_hash`), device memory
  management, event timing. Compiled only with `RF_HAVE_CUDA=1`.
- `Makefile` (`test`, `lib`, `lib-cuda`, `bench`, `ptxas-info`) and
  `CMakeLists.txt` (`check_language(CUDA)` → optional, as `gpu/mlwe`).
- `gpu/rootfind/rs/` — a **separate** cargo crate `gpu-rootfind`
  (`[workspace]` root of its own; depends on the root crate as
  `crypto_lib`). It holds the raw FFI, the safe wrappers `decompose64` /
  `decompose_p256`, the FFI tests, and the `gpu_rootfind_validate`
  harness. Its `build.rs` compiles the emulation library with the `cc`
  crate (no CUDA needed), or links a prebuilt `librootfind_cuda.a` +
  `cudart` when `GPU_ROOTFIND_CUDA_LIB_DIR` is set. The root crate
  deliberately carries **no** feature, dependency or build script for the
  oracle: cargo fingerprints a package with a build script by walking
  every file under it, and for the whole repository that walk fails on
  the CI runners (gix "Failed to update the excludes stack"); and an
  optional *path* dependency must exist at resolution time even with its
  feature off, which breaks the 14 workflows that build the root crate
  from a sparse checkout. A nested crate costs the root crate nothing.
- `src/cryptanalysis/gpu_rootfind.rs` (root crate, pure Rust):
  `Relation`; `cpu_reference64` (the same sweep on `Fp64::batch_inv` +
  `FastFactorBase::lookup`); `verify_relation64` (group re-addition);
  `synthetic_instance` (public known-answer fixture generator).

## 2. Build and run

### 2.1 This machine (no CUDA): what was run

```bash
make -C gpu/rootfind test                                  # C++ self-test: ALL PASS
cmake -S gpu/rootfind -B /tmp/build-rootfind && cmake --build /tmp/build-rootfind \
  && ctest --test-dir /tmp/build-rootfind                 # 1/1 passed (CUDA NOTFOUND → emulation only)
cargo build --release                                      # root crate: 0 build scripts, no C++
cargo test --release --lib gpu_rootfind                    # root: CPU reference test
cargo test --release --manifest-path gpu/rootfind/rs/Cargo.toml   # FFI + emulation vs reference: 4 passed
cargo run --release --manifest-path gpu/rootfind/rs/Cargo.toml --bin gpu_rootfind_validate   # PASS (§4)
cargo test --release --lib                                 # full suite (see PR)
```

### 2.2 RTX PRO 6000 Blackwell (sm_120, CUDA 12.8+): UNTESTED, exact commands

```bash
# 1. kernels + driver, native SASS for sm_120 plus PTX fallback
make -C gpu/rootfind bench ARCH=sm_120            # -gencode arch=compute_120,code=sm_120
                                                   # -gencode arch=compute_120,code=compute_120
make -C gpu/rootfind ptxas-info ARCH=sm_120       # registers / spills / smem per kernel

# 2. device vs host-emulation on identical inputs -- run this FIRST
cd gpu/rootfind
./bench selftest 61 1024 4096      # every (lookup, chunk) config; exit 0 iff all agree
./bench p256 256 1024              # P-256 instantiation
./bench throughput 61 65536 65536  # cells/s and targets/s; falsifies §3

# 3. the Rust harness against the device
make -C gpu/rootfind lib-cuda ARCH=sm_120
GPU_ROOTFIND_CUDA_LIB_DIR=$PWD/gpu/rootfind CUDA_PATH=/usr/local/cuda \
  BACKEND=cuda BITS=61 FB=4096 RAND_TARGETS=65536 \
  cargo run --release --manifest-path gpu/rootfind/rs/Cargo.toml --bin gpu_rootfind_validate

# alternative: CMake (CUDA auto-detected; archs default to 120-real;120-virtual)
cmake -S gpu/rootfind -B build-rootfind -DCMAKE_CUDA_ARCHITECTURES="120-real;120-virtual"
cmake --build build-rootfind && ctest --test-dir build-rootfind
```

Tuning knobs: `RF_BLOCK` (default 128), `RF_MIN_BLOCKS` (0 = ptxas
unconstrained), `chunk ∈ {16, 32, 64}`, `RF_SHARED_HASH_CAP`. Datacenter
Blackwell: `ARCH=sm_100`; Hopper: `ARCH=sm_90`.

## 3. Throughput model for sm_120 — ESTIMATE, not measured

Per `(target, i)` cell, 64-bit field, `CHUNK = 32`:

| Component | Field multiplications |
|---|---|
| batched inversion share: `3 + (62 sqr + 31 mul)/32` | ≈ 5.9 |
| two slopes `λ±` | 2 |
| two `x₃ = λ² − …` | 2 |
| `y₃` (hits only, ≈ `2N²/p` of cells) | ≈ 0 |
| **total** | **≈ 10 mul + ~8 add/sub + 2 probes** |

One Fp64 Montgomery multiplication on a 32-bit integer datapath is three
`IMAD.WIDE.U32` for the low product, four plus carries for `__umul64hi`,
the same again for `m·p`, plus adds and a select: ≈ 25 SASS instructions.
So a cell is ≈ 250 arithmetic instructions plus the probe (~15 for a
shared-memory hash hit/miss; ~16 dependent 8-byte loads for the sorted
table at `N = 2^16`) and loop overhead: call it **≈ 300 instructions/cell**.

Published RTX PRO 6000 Blackwell parameters: 188 SMs, 128 lanes/SM,
boost ≈ 2.4–2.6 GHz. Integer issue is the bound; taking 64 IMAD-class
lanes/SM/clk (the Ada/Hopper figure — Blackwell's unified INT32/FP32
cores may double this) gives

```
188 SM × 64 lanes × 2.4 GHz ≈ 2.9×10¹³ int-instr/s
÷ 300 instr/cell            ≈ 1×10¹¹ cells/s          (upper, perfect issue)
                              5×10¹⁰ cells/s          (with IMAD contention / 50 % issue efficiency)
```

→ **≈ 0.8–1.5×10⁶ targets/s at `|F| = 2^16`**, **≈ 5–10×10⁷ targets/s at
`|F| = 2^10`**. The emulation on this machine measured
1.0–1.9×10⁸ cells/s (multi-threaded, §4), so the estimate is a
300–1 000× ratio. P-256 on 64-bit limbs: ≈ 31 mul/cell at `CHUNK = 16`
(the 384-multiplication Fermat inversion dominates) × ~350 instructions
per SOS multiply ≈ 1.1×10⁴ instr/cell → **≈ 2–3×10⁹ cells/s**; the
32-bit CIOS field from `gpu/ecc` should improve that ~4×.

Occupancy (derived, to be checked with `ptxas -v`): Fp64 at `CHUNK = 32`
holds 64 u64 of inversion scratch = 128 registers + ~40 temporaries ≈
170 → 12 warps/SM of 48 on a 64K-register SM; `CHUNK = 16` roughly
doubles that at +30 % multiplications; `CHUNK = 64` spills. P-256 at
`CHUNK = 16` is already at the 255-register ceiling, hence its default.
Memory is not the bound: a cell reads 16 B of warp-uniform factor base
(L1) and nothing else on a miss.

## 4. Validation results (host emulation, this machine)

Harness: `gpu/rootfind/rs/src/bin/gpu_rootfind_validate.rs`. For each run: an `a = −3`
prime-order curve from `find_a3_curve(bits, 3)`, the small-x factor base
`FastFactorBase::build`, targets = positive controls `±F_i ± F_j` +
16 degenerate `±F_i` + random on-curve points + the `R = aG + bQ` behind
every 1-/2-entry relation `collect_relations` emitted on that curve
(200 000 targets scanned); then all 6 configurations of the oracle. Checks
per configuration: every relation re-adds in the group; relation **set**
equals `cpu_reference64` (the `prime_fast` sweep) on **every** target;
every control present; no duplicates. Apple-silicon host, emulation
backend, `BACKEND=auto` (resolved to `emulate`).

| Run | `p` bits | `|F|` | Targets (pos / degen / random / production) | Reference relations | Configurations | Result |
|---|---|---|---|---|---|---|
| `validate_40bit_fb4096` | 40 | 4 096 | 8 723 (512 / 16 / 8 192 / 3) | 531 | 6 | **PASS** 0/0/0/0/0 |
| `validate_56bit_fb16384` | 56 | 16 384 | 4 368 (256 / 16 / 4 096 / 0) | 272 | 6 | **PASS** 0/0/0/0/0 |
| `validate_61bit_fb4096` | 61 | 4 096 | 8 720 (512 / 16 / 8 192 / 0) | 528 | 6 | **PASS** 0/0/0/0/0 |
| default run (48-bit, `|F| = 1 024`, console only) | 48 | 1 024 | 2 320 (256 / 16 / 2 048 / 0) | 272 | 6 | **PASS** |
| `cargo test … gpu_rootfind` (40 + 61 bit, `|F| = 300`) | 40, 61 | 300 | 2 × 241 | — | 6 each | **PASS** |

(0/0/0/0/0 = re-add failures / duplicates / missing vs CPU / spurious vs
CPU / controls missing.) Total: **21 811 targets and 1 331 relations** in
the three evidence files (≈ 24 600 targets including the console run and
the unit tests), all configurations identical to the CPU oracle. Emulation
throughput 1.0–1.9×10⁸ cells/s vs the single-threaded Rust reference's
6–7×10⁷ cells/s. JSON + console logs:
`research/gpu_rootfind_decomposition_20261012/`.

`gpu/rootfind/test_cpu.cpp` (`make test`) additionally checks: Fp64
mul/add/sub/inv/batch-inv against `__int128` reference arithmetic over
4 primes × 20 000 samples incl. `p−1`; P-256 mul/add/inv/batch-inv against
a 512-bit schoolbook reference and the published `R mod p`, `R² mod p`,
`m'` constants; oracle positive controls on two synthetic curves with
sorted/hash × 16/32/64 agreement; P-256 positive controls; and the
fail-closed API paths (overflow detection, even `p`, coordinate ≥ `p`, bad
chunk, CUDA requested without CUDA → `RF_ERR_NO_CUDA`). `bench.cu`'s host
code was also compiled as plain C++ (warning-free) and run without a
device (exits 2, "no CUDA device").

**What is not validated:** anything below the launch line in
`rootfind_cuda.cu` — grid/block mapping, shared-memory staging, device
memory copies, event timing, `__umul64hi` codegen — and the Blackwell
flags themselves. `./bench selftest` (§2.2) is the check that closes
that gap; until it has run, the kernels should be regarded as
unverified on hardware.

## 5. Caveats and follow-ups

- `nvcc` absent here: `.cu` files were not compiled by any CUDA compiler.
  Likely first-run issues: `__launch_bounds__` placement with templates,
  the `RF_LAUNCH` macro's `<<<>>>` inside a `switch`, and static
  shared-memory size vs `RF_SHARED_HASH_CAP` on a given toolkit.
- P-256 on 64-bit limbs is a correctness mirror; swap in
  `gpu/ecc/fp256.cuh` (`Field<ModP,false>` with `FP_PTX=1`) for throughput.
- Global hash table for `|F| > 2048` (currently sorted + binary search)
  would cut probes from `log₂N` to ~2 loads.
- Degenerate `x_R = x_i` cells are patched to `d = 1` inside the batched
  inversion; this is also how `collect_relations` handles them.
- The oracle emits only `j ≥ i`; callers wanting both orientations
  should not symmetrise (they are the same relation).

## 6. Claim-hygiene boundary

This is a **decomposition subroutine** — the inner oracle of index
calculus, on **public synthetic known-answer factor bases** of small-x
points, with every returned relation re-verified in the group and the
whole result set cross-checked against the existing CPU implementation.
It changes the constant in front of the per-target cost of
2-decomposition; the number of targets that must be tried, the
relation-collection stage, and the linear algebra are untouched, so the
**asymptotic group-operation count of index calculus is unchanged**. No
key-recovery, asymptotic sub-rho, or deployed-curve-security claim is made
or implied; the P-256 instantiation is a field-arithmetic mirror used for
validation against the crate's own P-256 code, not a statement about
P-256. Hardware figures in §3 are estimates pending `./bench throughput`
on an sm_120 device.

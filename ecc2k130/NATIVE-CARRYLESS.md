# Native carryless multiplication

`PACKED_CLMAD=1` selects native carryless multiplication for the packed CUDA
backend. Its general default is zero. The RTX PRO 6000 preset enables it
with CUDA 13.3.1, batch 32, 256 threads per block, minimum two blocks and
256-worker state tiles.

The device implementation uses `clmad.lo.u64` and `clmad.hi.u64` with zero
addends. Both halves of each raw 64-bit product feed the existing three-leaf
Karatsuba reconstruction. The high-three-bit terms, field reduction, walk
formula, scalar counts, reports and checkpoint representation are retained.
Native multiplication takes precedence over `PACKED_GENERATED_PRODUCT`;
that generated software implementation is selected when CLMAD is disabled.
Host compilation retains software arithmetic.

NVIDIA introduced the operation in PTX 9.3 and documents support for targets
`sm_80` and later. This implementation requires CUDA 13.3 or newer, and its
native execution and performance were tested on `sm_120` using CUDA 13.3.73.
[PTX carryless multiply-add specification](https://docs.nvidia.com/cuda/parallel-thread-execution/index.html#integer-arithmetic-instructions-clmad).

## Complete-walk measurement

The [controlled comparison](benchmarks/clmad/comparison.json) ran both
implementations on one RTX PRO 6000 Blackwell Server Edition. It kept the
compiler, state layout, population and launch geometry fixed. Each timed
sample completed **201,863,462,912 scalar walk updates**.

| Workload | Software control median B/s | Native median B/s | Ratio-of-medians gain |
|---|---:|---:|---:|
| Complete walk benchmark | 7.110440 | **8.703518** | **22.4048%** |
| DP34 collection | 7.015447 | **8.534325** | **21.6505%** |

Native benchmark samples ranged from 8.701017 to 8.705273 B/s; collection
samples ranged from 8.533961 to 8.535006 B/s. Every paired measurement favored
the candidate. Two warmups and the initial screen are excluded from these
medians. All six collections contained the same sorted record multiset:
5,149 records, 164,768 bytes and zero drops.

Both GPU arithmetic suites, full client replay/restart tests, common-state
checks and bidirectional checkpoint-resume checks passed before timing. The
[independent retained-result audit](benchmarks/clmad/comparison-review.json)
verifies all 17 timed rows, exact counts, source and binary bindings, complete
native code and result statistics. The checkpoint comparison uses the same
TILE256 layout across implementations; its historical helper name still
contains “cross-layout.” Temporary corpus payloads were hashed during the
run; this artifact retains their sizes and multiset hashes.

The selected steady update path falls from 4,204.5625 to 2,233.0625 non-NOP
instruction visits per scalar update. The candidate uses 122 registers per
thread with zero stack, local or shared bytes. These static observations do
not predict an instruction-for-instruction speedup or establish a hardware
ceiling. The measured complete-walk rate remains below 15 B/s.

The separate [published-command audit](benchmarks/clmad/native-audit.json)
ran `make audit-rtx-pro6000` from the committed public source. It measured
**8.671278 B/s** benchmark median (8.671204–8.671598) and **8.507975 B/s**
collection median (8.506510–8.509265). All six samples completed the exact
scalar budget with the native-mode marker; each collection retained 5,149
records with zero drops. The [audit review](benchmarks/clmad/native-audit-review.json)
binds the source, binary, compiler, runtime flags and counts. This separate
allocation validates the command and does not estimate an additional gain.
It retains corpus counts and sizes; content equality comes from the paired
comparison above.

## Build and verify

For the fixed RTX preset:

```sh
make bench-rtx-pro6000
make audit-rtx-pro6000
```

For other existing build configurations, add `PACKED_CLMAD=1` to Make or
set `ECC_PACKED_CLMAD=1` for Modal, using CUDA 13.3 or newer. Other flags retain
their normal meanings. Invalid flag values and unsupported compiler/device
targets fail during validation or compilation.

The build cache, benchmark identity and audit record the requested mode.
Runtime output contains `packed native carryless multiply: 0|1`, with a
separate arithmetic-test marker. Missing, duplicated or mismatched markers
invalidate a result. A cached software build cannot satisfy a native request.

```sh
make test-clmad
make check-cli
```

The first target checks preprocessing guards and runs the independent host
field suite with CLMAD and generated-product selection enabled, exercising
the software fallback and selection precedence. Device correctness is
covered by the CUDA arithmetic and complete-client tests in the audit.

The [public-code calibration](benchmarks/clmad/public-code-identity.json)
binds the optional public implementation to the measured prototype: both
control and native cubins are byte-identical under the pinned compiler.

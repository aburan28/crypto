# Generated packed polynomial products

`PACKED_GENERATED_PRODUCT=1` selects a generated straight-line implementation
of the existing packed GF(2^131) polynomial product. Its default is zero,
retaining the handwritten product and reducer. Generated products require
`PACKED_DIRECT_REDUCE=1`. The RTX PRO 6000 preset enables generated products
and selects CUDA 13.3.1; other general defaults remain unchanged.

The generated function keeps the existing 144 widening-product events and
their coupled low/high outputs. It schedules native input preparation, the
top-three-bit correction terms and field reduction as one complete field
graph. Both outputs of the shared-operand polynomial helper use this function.
The normal-basis multiplier used inside inversion is unchanged. The walk
formula, batch grouping, scalar iteration count, distinguished-point rule,
report records and checkpoint format remain unchanged.

## Reproduction

From the `ecc2k130` directory:

```sh
python3 codegen/genpackedproduct.py --check
make test-generated-product
make test-packed-network
make bench-rtx-pro6000
make audit-rtx-pro6000
```

The generator uses only the Python standard library and reconstructs the
header deterministically. It does not load archived optimizer results or
private temporary files. `--check` verifies the checked-in header without
writing it. The generated function body matches the GPU-tested experimental
G1 schedule, apart from its public function name and leading comments.

The Modal environment variable is `ECC_PACKED_GENERATED_PRODUCT`, accepting
only `0` or `1`. Build reuse includes this selection. Packed measurements
require a unique runtime marker matching the requested mode; arithmetic
validation has a separate compiled-mode marker. Missing, duplicated or
mismatched identities invalidate the measurement.

## Controlled GPU result

The [raw comparison](benchmarks/generated-product/comparison.json) and
[independent review](benchmarks/generated-product/comparison-review.json)
cover one NVIDIA RTX PRO 6000 Blackwell Server Edition with 188 SMs and
driver 580.95.05. The control uses CUDA 13.0.48; the generated candidate uses
CUDA 13.3.73. Both use native sm_120 code, batch 32, 256 threads/block,
minimum two blocks/SM, 192,512 workers and 1,024 steps across 32 launches.
Each timed sample completes **201,863,462,912 scalar updates**.

| Workload | Native CUDA 13.0 median | Generated CUDA 13.3 median | Change |
|---|---:|---:|---:|
| Benchmark, B scalar updates/s | 6.860749 | **6.924275** | +0.925934% |
| DP34 collection, B scalar updates/s | 6.764850 | **6.819016** | +0.800698% |

Benchmark ranges across three repetitions are 6.858264–6.864944 B/s for
control and 6.923126–6.925822 B/s for the candidate. Collection ranges are
6.763079–6.765436 and 6.816272–6.820232 B/s. Every paired repetition favors
the candidate. Two warmups and three initial screening samples are retained
separately from these confirmation medians.

All six collection files contain 5,149 records, 164,768 bytes and zero drops.
Their sorted record hash is
`ab237b6352380547fd37fcdc9e2aa83f6ae1a718b842590b9a19e4224e5d336b`.
The source/binary/complete-code bindings, both arithmetic and client suites,
64-slot report/replay checks and 16,384-slot whole-state comparison passed.
Every one of the 17 timed samples completed its exact update budget without
a recorded error or timeout.

This is a modest improvement for the combined source, compiler and linked
runtime change. It does not establish an isolated generator speedup, a gain
on other GPUs/configurations, or the 15-billion-iterations/s target. Absolute
rates from other GPU allocations must not be used to infer another gain.

Raw artifact SHA256:
`87e658a3cdc685917c6451387735d09166f2a3c497bb611d7f160c6089ce26a4`.
Independent review SHA256:
`33f239ce5aa0c1fcf91d8893741239a8f64024d10749a6b82169dd4fb8309ea1`.

## Published-command audit and code binding

The [native preset audit](benchmarks/generated-product/native-audit.json)
ran `make audit-rtx-pro6000` from the public implementation at commit
`25db47d7ce6c0bc30ce21f60bb5186f9d2f9bf38` on a separate GPU allocation.
The [source manifest](benchmarks/generated-product/native-source-manifest.json)
retains all 109 input file hashes. Subsequent publication changes add only
documentation and evidence; the measured implementation remains unchanged.

| Workload | Median B scalar updates/s | Range across three repetitions |
|---|---:|---:|
| Complete walk benchmark | **6.960528** | 6.937184–6.963077 |
| DP34 collection | **6.809915** | 6.807981–6.817479 |

Every repetition completed 201,863,462,912 scalar updates with generated
product mode 1. Each collection recorded 5,149 points, 164,768 bytes and zero
drops. The GPU arithmetic suite and full client replay, restart, guard and
resume checks passed before timing. This audit retains collection counts
and sizes, but not collection content hashes; the paired comparison above
provides the matching-multiset evidence. The native audit verifies the
published command and is not an additional paired estimate of the gain.

The native source digest is
`d18667225ba454e86179756a6b036b0bfbb721107b2fe313bc4a6a1b1af49be4`;
the linked client binary hash is
`e1c2f918b0b0b4a28278346b3de91435e0bbbd6e34405383175eb03ae5aad76c`.
Raw native audit SHA256:
`c9686d2eaa20282acedbeef45ea0dbc1548ff2365fa7e7cbb11428b5ac221f72`.

A separate CPU compiler check compiled the public packed kernels with the
same CUDA 13.3.73 wrapper used for the tested G1 prototype. The complete cubin
and PTX are byte-identical, including all walk/init instructions and
encodings, helper placement and resource declarations. The public cubin hash
is `5b4d9ba9cdef7543563dada4aea2df4d2758ae21673a420aa85b2a369dfb88ac`.
The [raw compiler record](benchmarks/generated-product/public-code-resources.json)
and [independent review](benchmarks/generated-product/public-code-review.json)
retain commands, hashes and complete disassembly. This binds the public
device implementation to the measured prototype; it does not assert equality
of the complete linked host binaries or provide another throughput result.

## Compiler and correctness scope

Under CUDA 13.3.73, the same-source comparison reports 4,284.875 common-path
instruction visits per scalar update for handwritten products and 4,249.5
for the generated schedule. The ordinary and paired helper bodies decrease
from 564/1,100 to 556/1,087 instructions. The generated walk uses 128
registers, a 16-byte stack and no shared memory. Steady logical local traffic
is 0.5 bytes/update, plus 16 bytes/thread once per launch. These are compiler
and path-count observations, not cycle or DRAM-traffic estimates.

The independent derivation checks the complete 131-bit output map over
9,376 abstract input-bit columns, retaining carry-suppression masks and the
unused-high-bit convention. The generated code has been checked against all
17,161 canonical field-basis pairs and dense/edge cases, including host
Clang/UBSan checks of the actual header. Host kernel/reference fixtures cover
reports, guard restarts, zero denominators/numerators and exact normalized
state equality. GPU validation precedes the timed confirmation above.

The additional XOR-sharing experiment is not part of this feature. Its
smaller offline expression count increased live storage and did not improve
the compiled helper schedule. The checked-in generator emits only the
validated G1 schedule.

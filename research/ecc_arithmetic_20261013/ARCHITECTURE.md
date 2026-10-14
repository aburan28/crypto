# ECC arithmetic implementation map, 2026-10-13

**Status:** implementation analysis and test plan. The only new measured
kernel in this change is the binary-field CPU multiply in
[`src/binary_ecc/f2m.rs`](../../src/binary_ecc/f2m.rs). The backend options
below are candidates until built and measured. I treated the missing major
backend in the request as NVIDIA CUDA and included AMD HIP/ROCm as another
relevant platform.

## Separate workloads before tuning

| Boundary | Inputs and output | Measurement |
| --- | --- | --- |
| Field | Canonical elements, exact modulus; add/multiply/square/invert | latency and sustained throughput; include conversion and allocation where public API does |
| Point | Same curve model, projective/affine convention, identity and exceptional cases | add/double/normalize, with the coordinate conversion charged |
| Scalar | Same point, scalar distribution and secrecy policy | variable-base and fixed-base separately; cold precomputation and reusable precomputation separately |
| Batch | Exact number of independent points/scalars | wall latency, operations/s, bytes moved, setup and verification |
| ECC2K walk | Same curve, subgroup, partition, and distinguished-point rule | useful checked steps/s, not just field products or kernel iterations |

The current binary-field panel stops at the first row. A world-best result
cannot be inferred from a row at another boundary. For fixed-base tables,
report both cold-start cost and the exact number of uses that amortize the
table. Secret-scalar paths require constant-time control flow and memory
access; public rho walkers can use a different policy, which must be named.
The existing portable `clmul64_soft` path has an operand-dependent skip and
table index. This change does not establish constant-time behavior for
the binary-field API or any scalar operation built on it.

## CPU instruction paths

| Platform | Binary field | Prime field | Immediate implementation gate |
| --- | --- | --- | --- |
| ARM64 | PMULL/PMULL2 carryless 64-bit products; polynomial-basis sparse reduction; bit spreading for square | Scalar wide multiply, carry chains, curve-specific reduction; compare radix choices and interleave independent products | Runtime feature check, disassembly, ARM reference vectors, pinned P-core rounds |
| Intel x86-64 | PCLMULQDQ for one product, VPCLMULQDQ for lanes of independent products | BMI2 `MULX` with ADX `ADCX`/`ADOX` carry chains; curve-specific reduction | Feature-gated dispatch, equivalent fallback, ADX/BMI2 and VPCLMUL variants measured separately |
| Other CPU | Fixed-window software carryless multiply and standard integer limbs | Generic correct fixed-limb fallback | No assumption that a vector extension exists |

ARM documents PMULL as polynomial multiplication associated with FEAT_AES
([Arm architecture note, retrieved](https://developer.arm.com/community/arm-community-blogs/b/architectures-and-processors-blog/posts/webassembly-bitmask-operations)).
Intel lists PCLMUL and VPCLMUL intrinsics
([Intrinsics Guide, retrieved](https://www.intel.com/content/www/us/en/docs/intrinsics-guide/index.html))
and describes the `MULX`/`ADCX`/`ADOX` big-integer combination
([optimization manual, retrieved](https://cdrdv2-public.intel.com/821613/355308-Software-Optimization-Manual-048-Changes-Doc-2.pdf)).
Instruction availability alone predicts no runtime win: cross-limb XORs,
reduction, dispatch, register pressure, and clock behavior remain in the
measurement interval.

## GPU paths

| Platform | Candidate mapping | Principal risk and required measurement |
| --- | --- | --- |
| Apple Metal | 32-bit limb or bitsliced GF(2) kernels; one independent element/point per lane, with batch inversion only where it wins | Determine available integer width and thread execution width at runtime; measure command, copy, synchronization, and CPU reference verification, including small batches |
| NVIDIA CUDA | For supported devices, PTX `clmad.lo/hi.u64` for binary products; batched prime limb arithmetic and point updates | Check generated SASS and register occupancy; compare kernel-only to host-to-device plus device-to-host wall; validate every output |
| AMD HIP/ROCm | Test 32-bit limb and bitsliced carryless alternatives; specialize for the actual wave size | Query wave size rather than assume 32; measure occupancy, spills, and transfer costs |

Apple's [feature tables (retrieved)](https://developer.apple.com/metal/capabilities/)
show that Metal capabilities vary by GPU family, and its
[threadgroup guide (retrieved)](https://developer.apple.com/documentation/metal/creating-threads-and-threadgroups?changes=la)
requires obtaining `threadExecutionWidth` from the pipeline. NVIDIA's
[PTX ISA (retrieved)](https://docs.nvidia.com/cuda/parallel-thread-execution/)
documents `clmad` for SM 80 or newer in PTX ISA 9.3. AMD's
[HIP language guide (retrieved)](https://rocm.docs.amd.com/projects/HIP/en/docs-6.4.1/how-to/hip_cpp_language_extensions.html)
documents architecture-dependent wave size. GPU use is a batch-size decision;
a single scalar multiplication may lose to CPU after dispatch and transfer.

## Curve and representation decisions

| Curve/field | Baseline in this checkout | Candidate arithmetic, subject to verification |
| --- | --- | --- |
| ECC2K-130, GF(2^131) | Generic `F2mElement`; specialized walk kernels exist elsewhere in the repository | PMULL/PCLMUL three-limb Karatsuba, sparse reduction specialized for the exact modulus; square as bit spread; batch inversion; Koblitz Frobenius recoding for scalar multiplication |
| General GF(2^m) | Arbitrary-width `F2mElement` | Width-specific 2/3/4/6/9-limb multiplication crossover, CLMUL vs shift/XOR reduction, sparse/dense moduli kept distinct |
| NIST P-192/224/384/521 | General prime field `BigUint` path in the audited surface | Fixed-limb constant-time fields with modulus-specific reduction; P-521 merits a Mersenne fold; use verified generated arithmetic as a differential reference |
| NIST P-256 | Existing constant-time Montgomery field and projective point path | Compare fused multiply/reduce and point formulas against established implementations at identical secrecy and conversion boundaries |
| Curve25519 | `BigUint` field in `curve25519.rs` | Fixed radix 5x51 or 10x25/26 field, Montgomery ladder for X25519, RFC 7748 decoding and all-zero behavior |
| Ed25519 | Shares the `BigUint` field | Reuse verified fixed radix field, extended Edwards point formulas, separate signing and verification workloads and RFC 8032 encoding checks |

The field and curve identities follow the checkout's
[`curve_zoo.rs` (internal)](../../src/ecc/curve_zoo.rs) and
[`binary_ecc/curve.rs` (internal)](../../src/binary_ecc/curve.rs), and
[NIST SP 800-186 (retrieved)](https://csrc.nist.gov/pubs/sp/800/186/final).

### Prime reduction candidates

The NIST field primes have different reduction shapes, so a generic
Montgomery routine is a correctness baseline but not the sole candidate.
The exact modulus and candidate limb layout are separate choices:

| Field | Modulus identity | Candidate layout and reduction to test |
| --- | --- | --- |
| P-192 | `2^192 - 2^64 - 1` | 3 x 64-bit limbs; compare direct sparse folding with Montgomery |
| P-224 | `2^224 - 2^96 + 1` | 4 x 56-bit or 7 x 32-bit limbs; bound carry propagation before comparing to Montgomery |
| P-256 | `2^256 - 2^224 + 2^192 + 2^96 - 1` | Existing 4 x 64 Montgomery baseline; test fused reduction and verified generated output |
| P-384 | `2^384 - 2^128 - 2^96 + 2^32 - 1` | 6 x 64-bit limbs; compare direct folding and interleaved Montgomery carry chains |
| P-521 | `2^521 - 1` | 9 limbs with a short top limb or a radix chosen for carry bounds; Mersenne fold |
| 25519 | `2^255 - 19` | 5 x 51-bit on 64-bit CPUs and 10 alternating 25/26-bit on narrower SIMD/GPU lanes |

The local code does not currently provide these six fixed-limb paths as a
matched panel. For each candidate, prove or test bounds after multiplication,
carry and reduction, canonical output, and constant-time conditional
subtraction. A shorter reduction expression can still lose if it lengthens
the dependency chain. The cited [NIST parameter source](https://csrc.nist.gov/pubs/sp/800/186/final)
and [Fiat-Crypto examples](https://github.com/mit-plv/fiat-crypto/blob/master/Makefile.examples)
are retrieved references for these identities and candidate generators.

For GF(2^131), the exact local modulus is
`z^131 + z^8 + z^3 + z^2 + 1` in
[`f2m.rs` (internal)](../../src/binary_ecc/f2m.rs). PMULL/PCLMUL product
count is only one part of the cost; squaring can spread bits and fold the
result, while inversion should compare Itoh--Tsujii chains and batch
inversion at the actual batch size. A normal basis would make squaring a
cyclic shift, but conversion and multiplication costs must be charged before
considering it for this polynomial-basis workload.

The Curve25519 ladder and Ed25519 extended-coordinate behavior are specified
in [RFC 7748 (retrieved)](https://datatracker.ietf.org/doc/html/rfc7748)
and [RFC 8032 (retrieved)](https://datatracker.ietf.org/doc/html/rfc8032).
[Fiat-Crypto (retrieved)](https://github.com/mit-plv/fiat-crypto) provides
generated fixed-limb reference families for 25519, P-256, P-384, and P-521.

For binary Koblitz curves, Frobenius maps a point's coordinates to their
squares. A tau-adic scalar representation can trade some doublings for
Frobenius operations. It still pays for recoding, additions, field squares,
normalization, and precomputation, and its benefit depends on the exact curve
and workload. This is a scalar-multiplication path, separate from an
unknown-scalar ECDLP attack. A
[primary research treatment of endomorphism windowing (retrieved)](https://www.iacr.org/archive/eurocrypt2003/26560388/26560388.pdf)
and the checkout's Koblitz code are starting points, not evidence of this
implementation's speed.

## Reproducible comparison gate

1. Bind exact curve/subgroup/modulus, operation, input stream and public or
   secret scalar policy. Differential test all boundary values and RFC/NIST
   vectors; cross-check the final public API output, not only an internal
   product.
2. Build the same source revision and compiler options per arm; capture CPU
   model, microcode/OS, instruction flags, GPU driver, generated instructions,
   memory, frequency policy and isolation level. Include an A/A control and
   interleave A/B order.
3. Record every round, rejected sample, correctness failure and timeout.
   Apply the repository's measurement-quality gate before a wall comparison.
4. Compare against production libraries on the same host and boundary:
   [BoringSSL's benchmark binary (retrieved)](https://github.com/google/boringssl/blob/main/BUILDING.md)
   and [eBACS/SUPERCOP (retrieved)](https://bench.cr.yp.to/supercop.html)
   supply public baselines, but their published measurements from another host
   are not a matched local comparison.
5. For GPUs, separately report kernel time, allocation, transfer, launch,
   synchronization, verification, batch size and end-to-end latency. Include a
   same-host CPU arm and a no-op launch/transfer control.

The Mac panel failed step 3 because of host contention. After source-export
authorization, two Modal x86-64 field panels completed, but the workers
reported no CPU model, isolation, or A/A control and do not qualify a
wall-time comparison. There is no Metal, CUDA or ROCm result in this study.

# ECC arithmetic kernel study, 2026-10-13

The architecture and algorithm choices are detailed in
[ARCHITECTURE.md](ARCHITECTURE.md).
The unresolved measurement and implementation requirements are tracked in
[the problem record](../../problems/ecc-arithmetic-platform-measurement.md).

## Requested scope

The requested target is fast elliptic-curve arithmetic on ARM64, Intel CPUs,
Apple Metal, and a further major GPU backend (NVIDIA CUDA; AMD ROCm is also
relevant). Curves: NIST P-192, P-224, P-256, P-384, P-521; ECC2K-130;
Curve25519; Ed25519; general prime and binary fields; and Koblitz Frobenius
arithmetic. A worldwide-best claim requires a matched external comparison at
the same operation boundary and security/validation requirements.

## First measured slice

This study starts with the existing generic `F2mElement::mul` used by binary
curves. The candidate reduces 64-bit carryless products for two- and
three-limb elements. The three-limb case includes the degree-131 field used by
the ECC2K arithmetic code. On ARM64 it dispatches to PMULL, on Intel to
PCLMULQDQ, and otherwise to the existing portable path.

The small Cargo driver in this directory compiles the repository's actual
`src/binary_ecc/f2m.rs` and uses 256 deterministic input pairs per
field, a fixed seed `0x0123456789abcdef`, 200,000 products per round, five
rounds per field, and the same order for both builds. Before timing, each pair
is compared to `schoolbook_mul`. Retain every round as CSV, including any
failed run. The measured boundary is public `F2mElement::mul`, including its
allocation, dispatch, carryless product, and reduction. It does not measure
point addition, scalar multiplication, Frobenius chains, a DLP walk, or a GPU.

Build the baseline before editing `src/binary_ecc/f2m.rs`. Rebuild the same
driver after editing, then run `src/bin/compare.rs` with the baseline source
frozen under `baseline_src/` and the candidate imported from the live source.
It alternates AB and BA within one process and checks both against the same
bitwise reference. Run correctness tests on the candidate. A wall-time
comparison from a shared macOS host is descriptive; temperature, scheduling,
core placement, and other jobs are uncontrolled. Keep raw per-round data and
host/build metadata, and report regressions as well as gains. No benchmark row
is promoted to a general speedup claim from this panel alone.

## Work required beyond the first slice

| Track | Present implementation or measurement gap | Next gate |
| --- | --- | --- |
| ECC2K-130 | Generic binary field and specialized walk kernels have different accounting boundaries | Same-point, same-workload point/walk panel with reference answers and complete CPU/GPU transfer costs |
| NIST P-192/224/384 | General `BigUint` field/point paths remain; no matched optimized fixed-limb panel | Add fixed-limb field backends, differential vectors, and point/scalar benchmarks per curve |
| NIST P-521 | Nine-word field add/sub/mul/invert passed 1,000 local oracle rows; no point integration or timing | Connect to P-521 point path, validate scalar results, then run isolated same-boundary comparisons |
| NIST P-256 | Constant-time fixed-limb backend exists | Compare matched field and point operations to established implementations on each CPU |
| X25519 / Curve25519 | Five-limb field and ladder implemented; local differential checks passed, release-library and isolated timing gates open | Complete Cargo/RFC protocol suite, machine-code timing review, then isolated same-host A/A and A/B panels |
| Ed25519 | Signing, verification and point field arithmetic still use `BigUint` | Reuse validated fixed radix field in Edwards points, validate RFC 8032 and subgroup/encoding behavior, then benchmark sign and verify separately |
| Koblitz Frobenius | Tau/Frobenius implementations and generic arithmetic have different workloads | Validate scalar recoding and full scalar multiplication against a common reference on each field |
| Intel | No Intel host in the local run | Pinned x86-64 host with PCLMUL/VPCLMUL and BMI2/ADX feature manifest; paired rounds |
| Metal | No matched Metal arithmetic kernel or data-transfer accounting in this panel | Device checked batch kernel, CPU reference equivalence, dispatch and transfer timing |
| CUDA/ROCm | No matched GPU field/point panel in this study | Device checked batch kernel with occupancy, transfers, and end-to-end timing |

The X25519 row now has an implemented five-limb field and ladder in
`src/ecc/f25519.rs` and `src/ecc/x25519.rs`. Its
[local validation record](x25519_validation/README.md) retains 1,000 field
rows, 100 full-operation rows, two RFC vectors, and independent integer
checks. Release-library validation and an isolated paired timing panel are
still required. Automatic approval review rejected a newly frozen Modal
source export for this row before any remote launch; the
[exact bundle and denial](x25519_validation/MODAL_REVIEW_DENIAL.md) are
recorded separately from correctness.

The P-521 row now has an implemented nine-word Mersenne field in
`src/ecc/p521_field.rs`. Its [local validation record](p521_validation/README.md)
retains 1,000 seeded raw rows and independent Python checks of canonicalization,
addition, subtraction, multiplication and inversion. It has no point path or
timing claim. A [small source-bound Cargo harness](source_tests/README.md)
passed 13 focused tests across the 25519 and P-521 modules after the large
library test was stopped under host load.

The broad target remains open until every row has executable, verified
implementations and comparable measurements. New quantitative advantage claims
must pass the repository's claim and independent-review gates.

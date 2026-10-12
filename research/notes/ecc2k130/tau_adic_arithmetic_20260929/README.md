# Tau-adic arithmetic: useful scalar result, no established walk gain

2026-09-29 UTC. Completed bounded **CPU known-scalar arithmetic diagnostic**.
The m=83 and m=131 panels pass the preregistered arithmetic screening rule.
GPU throughput and end-to-end rho improvement remain **unmeasured**. No
production kernel, initialization, walk definition, or campaign changed.

## What the source audit established

Pinned source: [aburan28/crypto at 524d3877](https://github.com/aburan28/crypto/tree/524d38773a5ee5ba8f7ffe7e969241f2218dc1e4).
Exact inspected blob identities are in [source_audit.json](source_audit.json).

| Path | Existing operation | Consequence for this experiment |
|---|---|---|
| `ecc2k130/include/walk.h`, `ref.h` | One Frobenius-transformed point plus one point addition per ordinary non-table step | Already a sparse Frobenius/add expression; no binary scalar loop to replace with tau-NAF |
| `ecc2k130/include/tablewalk.h`, `packedtablewalk.cuh`, `packedkernels.cuh` | Select a signed Frobenius image of a precomputed table point, then add | Table selection and the addition remain; recoding a known scalar does not eliminate them |
| `ecc2k130/include/packedkernels.cuh`, initialization | Sum selected Frobenius images using 128 seed bits | Already uses a tau-polynomial representation; the integer-scalar diagnostic below does not measure this initialization |
| `gpu/ecc2k/koblitz.cuh` | Binary scalar multiplication for seeding/tests; source explicitly excludes it from the walk itself | A distinct known-scalar operation exists, but its fraction of runtime has not been measured |
| `ecc2k130/include/ref.h` | Binary scalar multiplication in reference/verification machinery | Potential arithmetic applicability; no measured end-to-end benefit |

This rejects **direct insertion of a long tau-NAF multiplication into either
audited hot loop** as unsupported by this evidence. It does not prove that
every future Frobenius-based walk or arithmetic optimization is impossible.

For rho generally, tau-adic scalar multiplication needs a suitable cheap
endomorphism. It is not a universal replacement for additions on arbitrary
curves. In particular the base-field Frobenius on an F_p-rational point is
the identity; this Koblitz result does not transfer directly to P-256 or
P-384. GLV/GLS are separate applicability questions.

## Measured arithmetic results

One unit: **microseconds per known-scalar multiplication**, median of seven
interleaved rounds, on the independent holdout (seed 20260930). Smaller is
better. All four reference implementations use identical affine formulas
over Sage field elements. Sage native multiplication is an additional oracle
and backend diagnostic, not a claim about the best available implementation.

| Variant | m=31 us | m=83 us | m=131 us | m=131 paired cost / binary NAF | Verified |
|---|---:|---:|---:|---:|---|
| Binary double-and-add | 159.22 | 659.40 | 1109.22 | 1.137 | Yes |
| Binary NAF control | 144.17 | 625.34 | 1004.01 | 1.000 | Yes |
| Unreduced tau-NAF | 127.74 | 478.77 | 776.67 | 0.776 | Yes |
| Reduced tau-NAF | 77.95 | 274.19 | 420.49 | 0.413 | Yes |
| Sage native backend | 501.93 | 1945.97 | 3517.25 | 3.561 | Yes |

The last numeric column is the median of *paired ratios*, not the ratio of
the two reported medians. Per-scalar conversion/reduction is timed. Curve,
field, and fixture construction is recorded separately and excluded from
this stage measurement. The annihilator is a shared curve constant.

No windows, native GPU backend, batch inversion, normal-basis GPU storage,
or fixed-base optimized scalar implementation is benchmarked here. This
Python/Sage experiment is evidence about its own arithmetic backend.

Deterministic arithmetic-expression counts explain the m=131 holdout:

| Variant | Additions | Doublings | Frobenius maps | Field products | Squaring expressions | Inversions | Recoding length |
|---|---:|---:|---:|---:|---:|---:|---:|
| Binary NAF | 43.29 | 127.92 | 0 | 342.42 | 299.13 | 171.21 | 128.92 |
| Unreduced tau-NAF | 85.17 | 0 | 255.67 | 170.33 | 596.50 | 85.17 | 256.67 |
| Reduced tau-NAF | 43.25 | 0 | 128.92 | 86.50 | 301.08 | 43.25 | 129.92 |

Counts are averages over the same 24 scalars. These columns have different
costs; they are **not summed into an instruction count or ECDLP S metric**.
Field additions, recoding integer work, allocation, and language overhead
are included in elapsed time but not in this selected operation vector.
Reduced tau-NAF removes about 75% of the inversions versus binary NAF in
this affine arithmetic. Unreduced recoding illustrates why reduction matters.

## Noise and preregistered decision

| Degree | Panel | Reduced tau / binary NAF paired cost | Largest A/A deviation from parity | Screen |
|---|---|---:|---:|---|
| 31 | Training | 0.546 | 17.32% | Pass |
| 31 | Holdout | 0.542 | 48.05% | **Fail: gain did not exceed observed noise** |
| 83 | Training | 0.509 | 46.79% | Pass |
| 83 | Holdout | 0.454 | 31.74% | Pass |
| 131 | Training | 0.394 | 19.22% | Pass |
| 131 | Holdout | 0.413 | 22.34% | Pass |

The primary m=83/m=131 screening criterion passes, but short batches on a
shared VM produced substantial timing variation. Do not treat these as
precise performance ratios or as a confidence-interval-qualified runtime
claim. The failed m=31 noise gate is retained. The run stopped after the
frozen workload; no favorable rerun was selected.

Host: Linux x86-64, Intel Xeon Platinum 8370C, family 6/model 106, KVM,
SageMath 10.9. Pinned to logical CPU 0. NUMA memory binding to node 0
succeeded and was queried back (MPOL_BIND, mask 1). One visible NUMA node;
CPU pinning does not reserve a core or reveal physical-host contention.
DIMM/DDR type is not exposed and remains unknown. No CUDA compiler,
`nvidia-smi`, or `/dev/nvidia*` device was available. Raw receipts contain
CPU features, cgroup limits, process names, and before/after load snapshots.

## Correctness and evidence

- 2,310 exhaustive point/scalar cases across E_0 and E_1 over GF(2^5),
  including infinity and torsion: 9,240 candidate/oracle equalities.
- 135 larger-field edge cases and 144 train/holdout cases: 1,116 more
  candidate/oracle equalities. **10,356 checked equalities total**.
- Independently verify every candidate against Sage point multiplication;
  check tau's quadratic identity, tau^m=1, exact digit reconstruction,
  nonadjacency, irreducible moduli, and subgroup membership of fixtures.
- All timed and counted reruns produce the same output digest. Six panels
  retain exact point/scalar inputs, seeds, hashes, every sample, and counters.
- m=31's order divided by four is composite, recorded explicitly. This is
  a scalar-arithmetic test, not a prime-subgroup rho benchmark. The m=83 and
  m=131 quarter orders are prime.
- [run-01.json](results/run-01.json) retains the first failed attempt: its
  exhaustive correctness tests passed, then NTL field-element serialization
  failed before measured panels. The single portability fix was to serialize
  polynomial coefficients; it did not change arithmetic or inputs.
- [run-02.json](results/run-02.json) is the completed run. The receipt verifier
  checks its current source/protocol hashes, ratios, counters, output digests,
  all six panels, the retained failed receipt, and null end-to-end fields.

The protocol and initial source were locally committed as `c28d3b7` before
execution, and preserved in remote commit
`49c87f2b46ca5f537f1edfc18635ca203b137fd2`. The serializer fix and failed
receipt were locally committed as `168d5ae` before the second run.

## Reproduction

Run from this directory. A new result path is mandatory; the frozen protocol's
original `run-01.json` destination is already occupied by preserved evidence.

```bash
python3 verify_receipt.py
sage -python benchmark.py --output results/local-reproduction.json
```

The implementation is variable-time research code. Correctness against Sage
does not establish suitability for secret-key operations in production.

## Decision and remaining questions

- [x] Audit both walk families and initialization at a fixed source revision.
- [x] Freeze the protocol, implement scalar controls, verify, and run both
  training and independent holdout panels at m=31/83/131.
- [x] Preserve the failed run, operation vectors, hardware/NUMA evidence,
  raw paired timings, source hashes, and explicit negative noise result.
- [x] Reject a claim that this scalar result already accelerates the hot loop.
- [ ] GPU known-scalar arithmetic comparison: not run here.
- [ ] Measured runtime fraction of eligible known-scalar work: unknown.
- [ ] Whole-walk or end-to-end rho benefit: not established.

For illustration only, if eligible scalar work were 1% of total time and
its candidate/control cost ratio were 0.413, Amdahl's law predicts only
about 0.59% total time saved. Both the runtime fraction and GPU ratio are
unknown, so the actual full-application fields remain null.

Classification: promising **arithmetic engineering diagnostic**, no
cryptanalytic advance and no end-to-end performance claim. This is not an
index-calculus experiment; no IC scoreboard numbers are changed.

Literature and the preregistered rules are in [PROTOCOL.md](PROTOCOL.md).

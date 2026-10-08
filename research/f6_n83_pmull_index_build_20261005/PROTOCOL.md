# n=83 F6 packed pair-index construction experiment

Registered before implementation and timing on 2026-10-05. The rejected
shared-field experiment in `research/f2m83_pmull_20261005` accelerated
full index construction but missed a small-query regression gate. This
candidate confines packed ARM64 PMULL arithmetic to ordinary additions
while building `F6SignedPairIndex`. It computes both coordinates with
the exact degree-83 polynomial `z^83 + z^45 + z^2 + z + 1`, keeps one
reference inversion per batch, and calls the reference group law on
exceptional cases. All other curves and CPUs use the existing builder.
The target-query method, lookup, sign check, and witness verification
are unchanged. A package-internal raw-word constructor may avoid
BigUint conversions when materializing the exact affine output.

Baseline: retained query branch head `8b6e6e395`, with the tested
binary and source hashes in `research/f6_n83_pmull_xonly_20261005`.
Freeze K0 curve `icv1-f2m83-tm6151469093347-debefd74`, standard
cofactor-projected bases at dimensions 8, 10, and 12 (258, 1,048,
and 4,054 usable points), and public T001
(x=`355fb5df7a905f16921eb`, y=`5900a390f42d290f1bbe`). The full
index has 8,219,485 unordered pairs and 4,108,723 signed-sum
representatives. Build matched baseline and candidate native release
probes from exact source snapshots with the same compiler and flags.

Before timing, compare packed full-point additions to the reference
group law on ordinary and exceptional n=83 inputs, including identity,
same-x, inverses, and repeated points. Compare small-base four-sum
search to exhaustive group enumeration. Run a full-base planted K0
witness control and replay all returned points independently. Preserve
fallback correctness on nonmatching fields and unsupported hardware.

Run three repetitions per arm at each small base, then full-base probes
in baseline, candidate, candidate, baseline order, with no builds
between arms. Preserve exit codes, stdout/stderr, source and binary
hashes, index-build and query costs, RSS, pair counts, and witness
replay. Retain this builder only if all exactness controls pass,
full-base index-build median falls by at least 20%, full-base query
median does not regress by more than 10%, and peak RSS rises by at
most 10%. Record small-base query costs as diagnostics; their code path
is unchanged and they are not a retention gate. Otherwise revert
runtime code and preserve a rejected patch and all outcomes.

The host is unisolated; CPU ratios remain exploratory stage diagnostics.
Index build is target-independent and excluded from one-target online
IC time. The K0 four-summand base's uniform-target relation coverage
ceiling is `4.662e-12`. This experiment cannot establish a complete
F6/F4/F5 or one-target IC speedup.

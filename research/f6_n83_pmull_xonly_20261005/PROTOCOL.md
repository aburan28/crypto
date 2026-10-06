# n=83 F6 packed PMULL residual experiment

Registered before implementation and timing on 2026-10-05. The retained
sign-folded x-only query spends its main work calculating millions of
residual x coordinates in the fixed 83-bit field. The candidate computes
ordinary batch residuals as packed `u128` field elements using ARM64
PMULL multiplication and exact reduction modulo
`z^83 + z^45 + z^2 + z + 1`. It keeps one reference field inversion per
batch, uses the reference group law for exceptional points, and retains
the existing full hash, sign, and group-witness checks. The new method
falls back to `solve4_xonly` unless the field and ARM64 AES/PMULL
capability are confirmed. It does not alter index construction.

Baseline: x-only branch head
`f30cc085ed0c0a947b9c501c438f644634034239`, with pair-index source
SHA-256 `3f2ef8bd0f2bb07dde0bdb7a9ad7f3da76e0ff1924a9ec04dbd200102ae4125f`.
Freeze the K0 n=83 curve `icv1-f2m83-tm6151469093347-debefd74`,
standard polynomial-subspace cofactor-projected bases at dimensions
8, 10, and 12 (258, 1,048, and 4,054 usable points), and public T001
(x=`355fb5df7a905f16921eb`, y=`5900a390f42d290f1bbe`). The full
index has 8,219,485 unordered pairs and 4,108,723 signed-sum
representatives. Build baseline and candidate native release probes
with the same compiler and flags from their exact source snapshots.

Before timing, compare the packed-field multiplication and squaring to
the reference field on deterministic edge values and at least 10,000
deterministic pseudorandom pairs. Compare packed residual x keys to
the reference batch on ordinary and exceptional n=83 points, and
compare four-sum search against exhaustive small-base enumeration,
including infinity, sign choices, repeats, witnesses, and misses.
On non-ARM or a CPU without the feature, check the fallback result.

Run one three-repetition small-base probe per arm, then full-base
probes in baseline, candidate, candidate, baseline order, with no builds
between arms. Preserve all exit codes, stdout/stderr, source and binary
hashes, build and query costs, RSS, pair counts, and witness replay.
Retain the packed path only if all exactness controls pass, full-base
query median falls by at least 20%, neither small-base query median
regresses by more than 10%, full-base build cost does not exceed baseline
by more than 10%, and peak RSS rises by at most 10%. Otherwise revert
runtime code and preserve the rejected patch and all outcomes.

The host is unisolated; all ratios remain exploratory stage diagnostics.
The four-summand base's uniform-target relation coverage ceiling is
`4.662e-12`. A faster exact miss does not establish faster complete F6,
F4/F5, one-target IC, or Pollard rho.

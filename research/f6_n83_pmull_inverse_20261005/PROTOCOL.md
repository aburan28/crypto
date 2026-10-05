# n=83 F6 packed inversion experiment

Registered before changing source or measuring on 2026-10-05. The
retained packed query and pair-index builder still invoke the portable
`F2mElement::flt_inverse` once per fixed-point row. This experiment
replaces only those inversions, inside the exact degree-83 PMULL path,
with the same Itoh–Tsujii exponentiation using packed multiplication
and squaring. The reference field, other curves, and CPUs keep their
existing behavior.

Baseline is #1399 head `adb3d8b6877419bfa14186d5ea865b5eb98d0a11`.
Its F6 and field source SHA-256 values are
`887677cbc6c91bd996caf92e79fb575d0af6295fbb1e4fb0eb4c9a7cefb2e619`
and `16ea7db9ff4792f54e3f1b21454a7a2d109b9a2e27245f93e27c1ec530b0559e`.
Freeze the K0 curve `icv1-f2m83-tm6151469093347-debefd74`,
cofactor-projected standard subspace bases at dimensions 8, 10, and 12
(258, 1,048, and 4,054 usable points), and public T001
(x=`355fb5df7a905f16921eb`, y=`5900a390f42d290f1bbe`).
The dimension-12 index has 8,219,485 unordered pairs. Keep the
`solve4_xonly_pmull83` query and signed-pair lookup otherwise unchanged.

Before timing, compare packed inverse against the reference field
inverse on zero, one, and 10,000 deterministic nonzero field elements;
verify `a * inverse(a) = 1`. Repeat exact geometry tests for ordinary
and exceptional pair additions, the exhaustive small-base four-sum
test, and a full-base planted K0 witness with independent group replay.
Preserve a nonmatching-field/unsupported-hardware fallback check.

Build matched native Rust release probes from pinned snapshots with
the same compiler and flags. Run three repetitions per arm on each
small base, then full-base baseline, candidate, candidate, baseline
with no builds between full arms. Preserve exit codes, stdout/stderr,
source and binary hashes, build/query times, RSS, pair counts, and
witness replay. Retain the candidate only if every exactness control
passes, full-base **query** median falls by at least 15%, full-base
index-build median rises by at most 10%, and peak RSS rises by at most
10%. Small-base costs are diagnostics, not gates. Otherwise revert
runtime code and preserve a rejected patch and all outcomes.

The host is unisolated; timing ratios remain exploratory stage
diagnostics. Index construction is target-independent and excluded
from the primary one-target online IC interval. Four-summand coverage
on this base is at most `4.662e-12` for a uniform target. This
experiment cannot establish a complete F6/F4/F5 or one-target IC/rho
speedup.

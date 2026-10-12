# Shared degree-83 PMULL field arithmetic experiment

Registered before implementation and timing on 2026-10-05. The retained
F6 packed residual query still builds its 4.1-million-representative
index with the general binary-field arithmetic and uses the general
Itoh-Tsujii inversion once per 256 residuals. The candidate adds an
ARM64 AES/PMULL fast path to `F2mElement::mul`, `square`, and
`square_k_times` for the exact irreducible polynomial
`z^83 + z^45 + z^2 + z + 1`. It uses packed `u128` operations and
fixed sparse reduction; every other field and CPU retains the existing
portable implementation. The F6 query algorithm, index, and witness
checks remain unchanged.

Baseline: packed-query branch head `8b6e6e395` and its tested
source/benchmarks in `research/f6_n83_pmull_xonly_20261005`. Freeze the
K0 curve `icv1-f2m83-tm6151469093347-debefd74`, the standard
cofactor-projected bases at dimensions 8, 10, and 12 (258, 1,048,
and 4,054 usable points), and public T001
(x=`355fb5df7a905f16921eb`, y=`5900a390f42d290f1bbe`). Build
matched baseline and candidate release probes with the same compiler,
flags, and exact source snapshots. The full index has 8,219,485
unordered pairs and 4,108,723 signed-sum representatives.

Before timing, compare candidate multiplication and squaring against
the independent schoolbook method on field edge values and at least
10,000 deterministic pseudorandom pairs. Check `square_k_times` and
inversion round trips. Run the registered K0/K1 curve and subgroup
controls and the F6 exhaustive small-base and planted full-base
witness checks. A portable fallback or requested unsupported hardware
must remain explicit. A failure, timeout, or OOM is preserved.

Run three repetitions per arm at each small base, then full-base probes
in baseline, candidate, candidate, baseline order, with no builds
between arms. Preserve all exit codes, stdout/stderr, source/binary
hashes, index-build and query costs, RSS, pair counts, and witness
replay. Retain the shared field path only if all exactness controls
pass, full-base index-build median falls by at least 20%, full-base
query median does not regress by more than 10%, neither small-base
query median regresses by more than 10%, and peak RSS rises by at most
10%. Otherwise revert runtime code and preserve a complete rejected
patch and all outcomes.

The host is unisolated; all CPU ratios are exploratory stage diagnostics.
Index build is target-independent and excluded from the primary one-target
online IC interval. The K0 four-summand base's uniform-target relation
coverage ceiling is `4.662e-12`. This experiment cannot establish a
complete F6/F4/F5 or end-to-end one-target IC speedup.

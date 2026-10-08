# n=83 F6 exact x bitmap prefilter experiment

Registered before implementation and timing on 2026-10-05. The retained
sign-folded x-only query still probes a large hash table for each residual
x key, generally missing on the public T001 target. The candidate adds a
2^25-bit (4 MiB) bitmap over the low 25 bits of each stored affine x key.
It tests this bitmap before the existing hash lookup. Collisions may cause
extra lookups but cannot omit a stored x key. Infinity bypasses the bitmap.
The full hash, sign check, and curve-group witness replay remain authoritative.
The existing full-coordinate `solve4` is unaffected.

Baseline: `codex/f6-n83-xonly-query-20261005` head
`e40451d5fafcf5acabc02e98a339aa9081d2711d`. Freeze the K0 n=83
curve `icv1-f2m83-tm6151469093347-debefd74`, standard polynomial-subspace
cofactor-projected bases at dimensions 8, 10, and 12 (258, 1,048, and
4,054 usable points), and public T001 (x=`355fb5df7a905f16921eb`,
y=`5900a390f42d290f1bbe`). The full index has 8,219,485 unordered
pairs and 4,108,723 signed-sum representatives. Build baseline and
candidate native release probes with the same compiler and flags. Check
every stored x key against the filter, run the exact x-key/reference
group-law controls, and compare four-sum search to exhaustive small-base
enumeration before timing.

Run one three-repetition small-base probe per arm, then full-base probes
in baseline, candidate, candidate, baseline order, with no builds between
arms. Preserve all exit codes, stdout/stderr, source and binary hashes,
index-build cost, query cost, peak RSS, pair counts, and witness replay.
Retain the prefilter only if all exactness controls pass, the full-base
query median falls by at least 10%, neither small-base query median
regresses by more than 10%, full-base build cost does not exceed baseline
by more than 10%, and peak RSS rises by at most 10%. Otherwise revert
the code and preserve the rejected patch and all outcomes.

This unisolated-host ratio is an exploratory stage diagnostic only. The
four-summand base's uniform-target relation coverage ceiling is
`4.662e-12`. A faster exact miss does not establish faster complete F6,
F4/F5, one-target IC, or Pollard rho.

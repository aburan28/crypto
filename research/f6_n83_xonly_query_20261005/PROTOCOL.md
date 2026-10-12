# n=83 F6 sign-folded x-only residual query experiment

Registered before execution on 2026-10-05. The retained sign-folded
four-summand index looks up a residual by its x-coordinate, then checks its
sign. For an ordinary no-witness query, calculating every residual
y-coordinate may be avoidable work. The candidate will batch-compute both
residual x-coordinates for each stored signed pair sum, probe the existing
x-key lookup, and calculate a full residual point only when that x key is
present. Exceptional same-x and infinity cases retain the reference group
law. Every returned witness must replay in the curve group.

The baseline is the sign-folded implementation at
`d9f11431a81315fa2612982fd930b405c73e30d5`. Freeze the K0 curve
`icv1-f2m83-tm6151469093347-debefd74`, the standard polynomial-subspace
cofactor-projected bases at dimensions 8, 10, and 12 (258, 1,048, and
4,054 usable points), and public T001 (x=`355fb5df7a905f16921eb`,
y=`5900a390f42d290f1bbe`). The full index has 8,219,485 unordered
pairs and 4,108,723 signed-sum representatives. Build both native release
probes with the same compiler and flags. Prove candidate x keys against full
reference group additions, including exceptions, and compare candidate
four-sum search against exhaustive small-base enumeration before timing.

Run one small-base three-repetition probe per arm, then full-base probes
in baseline, candidate, candidate, baseline order. Preserve all exit codes,
stdout/stderr, source and binary hashes, index-build cost, query cost, peak
RSS, pair counts, and witness replay. Retain the x-only method only if
all exactness controls pass, the full-base query median falls by at least
15%, neither small-base query median regresses by more than 10%, full-base
build cost does not exceed baseline by more than 10%, and peak RSS rises
by at most 10%. Otherwise revert code and preserve a complete rejected
patch and all outcomes.

This unisolated-host ratio is an exploratory stage diagnostic only. The
four-summand base's uniform-target relation coverage ceiling is
`4.662e-12`. A faster exact miss does not establish faster complete F6,
F4/F5, one-target IC, or Pollard rho.

# Volcano IC: controls and 99 vertices across all 51 strata

Source: [native execution study](../../research/full_volcano_graphwide_ic_20261017/execution/README.md).
The dev pilot has 48 attempts (30 verified, 14 exhausted, 4 timeouts); its bounded
survey checkpoint now covers 99 vertices and all 51 strata, with 1,584 attempts
(1,065 verified, 422 exhausted, 97 timeouts). Every successful/exhausted
deterministic report has a bound replay receipt. This is not full
graph completion or admitted optimized timing. Original A--E and finite-field
IC requirements remain open.

The historical reports and SVG/150-dpi PNG graphs are copied here from
the study's write-once run reports. Native CSV/JSON statistics retain every
phase, status, stratum and present-value denominator; no missing timeout counter
is treated as zero. No best-method/frontier or scaling claim is promoted.

The new [all-strata report](all-strata-99-report.pdf) and
[SVG](all-strata-99-yield.svg)/[150-dpi PNG](all-strata-99-yield.png) show the
native per-stratum min/max/mean and both yield denominators. The checkpoint has
57,176 relations from 2,223,204 source trials and 2,075,145 actual oracle calls;
the old reports' "decompositions" label meant source trials, including queries
skipped before calling the oracle. Frozen old reports are retained unchanged.

The separate release control/pilot have 56 identical completed deterministic
comparisons against dev and four censored pairs. These are parity checks, not
an admitted runtime speedup. External comparison: not run; no fastest-method
claim is made. The optimized whole-queue continuation uses its own campaign.

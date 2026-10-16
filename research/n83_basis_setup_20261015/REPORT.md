# Degree-83 IC basis/setup: final publication of the attempted bounded study

C3 reduced actual cold process mean from K600's1223.7821487055s to
611.139215194875s on eight fresh matched targets. Optimized rho averaged
2.10772421875s. Native C3/rho is289.952171995877, conditional95% interval
[245.967289930960,375.396780301831]; native C3/C0 is0.499385626634
[0.484205609743,0.514186089767]. All IC/rho ratios below1 favor IC. Thus this
measured cold boundary does not produce an IC-over-rho end-to-end win.

[22-page inspected addendum](reports/final-20261016/n83_basis_setup_addendum.pdf) ·
[109-page comprehensive PDF](reports/final-20261016/unified_curve_basis_report_basis_setup.pdf) ·
[native final binding/profiles/statistics](reports/final-20261016/native_final.json) ·
[cold native statistics](reports/run-20261015T1842Z/native_statistics.json) ·
[cohort stage report](reports/run-20261015T2311Z/REPORT.md).

## Actual coverage and failures

| Stage | Attempted / completed / failed processes | Canonical answers | Full rank/log audits | Identical rho A/A work pairs |
|---|---:|---:|---:|---:|
| Training |18 /16 /2 |16 |8 |3 |
| Fresh cold |26 /26 /0 |26 |16 |1 |
| Reused8-target cohorts and final A/A |11 /8 /3 |50 |6 |1 |
| Actual total |55 /50 /5 |92 (72IC,20rho) |30 |5 |
| Prospective plan |55 /55 /0 |118 (74IC,44rho) |32 |5 |

Two C4 training processes failed a12884901888-byte allocation under the common
12GiB address-space policy. Their attempted costs and partial traces remain;
C4 was ineligible for selection. C3 was selected using both complete training
runs (mean622.6626231075s) and committed19664f02f before holdout generation.

All three cohort rho processes returned2 without answers. The frozen rho
worker accepts one or32 public points, not the predeclared eight. This was an
implementation-contract error that should have been detected before timing,
not evidence of expensive rho search. Failed launches are never rho solve
costs. The sealed v1 checker preserved its failed receipt; the full-stage
analyzer correctly refused to produce complete cohort statistics. No solver,
training, selection, generator or prior audit was rerun; no partial directory
was retried or downscaled. [Named unresolved gap](../../problems/n83-eight-target-rho-contract.md).

The six verified cohort IC processes have means C0=1150.628969694667s and
C3=591.192669254667s; native C3/C0=0.513799569475
[0.505576379485,0.522482172825]. This uses20,000 paired draws over only three
preparation clusters (seed202610157802), and is an IC-basis comparison only.
Cohort IC/rho, rho online rows and their intervals remain unavailable/null.
No historical32-target or cold rho times fill the missing cohort cells.

## Boundary, bases and accounting

Exact `icv1-f2m83-t6151469093347-9eba6cd5`, a=b=1,
GF(2^83): z^83+z^7+z^4+z²+1 (tail149), r8569786107849059,
G=(6614244651398184860292015,8868260410029249391539529), lambda2998004823950457.
No tail45/a=0 substitution. Synthetic public targets; validation scalars remain
outside measured workers. Common CPU3/oneOSthread/12GiB/2100s envelope.
Rho Fold/reference-u64/256lanes/DP6/escape-off/capfactor2000. Its seed-specific
G-only jump preparation and target-dependent DP/walk state are charged afresh;
no cross-target collision sharing.

C0 is retained K600 wide; C1 identical K600 dense; C2 next600 public-scan signed
orbits disjoint from K600 (construct1200 then drop/reindex600, all charged);
C3 public-scan prefix1271 dense; C4 prefix1798 dense. Dense table packing is an
engineering representation, not a different mathematical base. The first600
draft was an exact K600 duplicate and was corrected before timing. All bases
reconstruct relation logs; no supplied known logs are substituted.

Training targets exclude1911 old/base signed orbits; fresh holdouts generated
after committed selection exclude1913, including both training targets.
Distinct rank seeds and prospective orders are in immutable protocol/matrices.
Cold primary20k paired target/preparation-block bootstrap seed202610157801;
conventional midpoint medians/type7 p90/p95, mean costs/CPU/peakRSS and all
attempts appear in the addendum and native JSON. All eight cohort target rows
show all three available IC repeats and explicitly missing rho repeats.

Aligned online receipts in both methods include generic Q validation,
solve/recovery, generic scalar replay and answer serialization/write/flush;
exclude input-file parsing and own receipt emission. Cold secondary native
mean aligned query ratio C3/rho=0.1155785234636542 (full precision in native_final.json).
Actual process costs remain primary; query clocks do not remove setup.

## Preparation profile and independent verification

Thirty full transcripts independently rebuild every row/rank gain/backsolve,
verify each representative[log]G, anchor/four-point witness/intermediate and final
scalar/log reconstruction. Checker wall/CPU/RSS is separately wait4-metered,
excluded from solver times and retained in each stage audit directory.
Native final aggregation checks immutable raw/execution/matrix/input/audit
bindings and recursively compares all five A/A outcomes. Raw SHA tamper,
recursive A/A mutation and planned-versus-actual coverage regression guards pass.
The full-stage v1 source/binary seals still pass unchanged.

Index build and relation search dominate the recorded profiles; elimination
and log backsolve are small. Native JSON exports every completed-IC attempt's
probes, states, useful relation, rank before/after, rank gain, exclusive phase
times and endpoint RSS/HWM. C3 uses a4GiB packed_u128 root table (16-byte slots)
and6GiB wide_u128_roots states (48-byte slots). C4 failed under this particular
wide-state allocation policy; packed-state/larger-base impossibility is not
established. Endpoint samples are not exclusive phase peaks; producer counters
and timings are not independently re-executed search timing measurements.

## Inference, canonical views and preservation

One guest, only three preparation clusters; physical-host exclusivity and
independent-host replication are unestablished. PSI unavailable and unmovable
daemons recorded. Memory-limit-hit remains null. Original retained-K600
selection/construction, OS-cache flushing and instrumentation calibration remain
unpriced. No S/L2/global-optimum/exponent/asymptotic/deployed-security claim,
1685-target experiment or projected winner. Never multiply gains across panels.

The scoreboard/progress and generated leaderboard receive local diagnostics,
with historical results retained. Boundary ledger/admitted-best/floor rankings
stay unchanged: no calibrated S, sealed ecbench L2 or independent-host admission;
cohort reference is missing. Canonical graph check details are in checks/final.

The old87-page PDF remains byte-identical. Its text and all pixels at72dpi
(Fitz) and independently144dpi (Poppler) match the first87 pages of the updated
109-page comprehensive PDF. All22 new pages were read; no empty/out-of-bounds
pages. The standalone addendum matches the appended pages at both resolutions.
[Page-by-page preservation receipt](reports/final-20261016/PDF_PRESERVATION.json).
Prior70-page PDF/139-launch rho evidence and old reports remain historical.
Python only formats/presents native values, renders/inspects and packages files.

## Requirement status and review

| Requirement | Status |
|---|---|
| Relation/rank/log profile and transcript export |Verified for all30 completed IC processes; failure partials preserved |
| Five-candidate training/frozen fresh holdouts |Verified; C4 ineligible failures retained |
| Matched full-cost cold comparison |Verified, rho cheaper |
| Matched8-target batch comparison |Partial/unresolved: three rho processes rejected input |
| All attempted evidence/limits published |Published with no unqualified incomplete-panel winner |
| Inspected final PDFs/graphs/old87 preservation |Verified, source and receipts linked above |
| Current applicable CI and required human review |Merge gates; no merge or reviewer request |

PR1899 remains dependent on unmerged parent1828. Required human review including
Cursor remains outstanding; no reviewer is requested without owner naming one.
Publication records the failed cohort obligation rather than silently marking
its experimental acceptance complete. A correction requires a new version and
prospective study, never editing these frozen workers/records.

## Evidence delivery and final validation

[Evidence ZIP](delivery/n83_basis_setup_evidence.zip),
[exhaustive hash/reconstruction manifest](delivery/HASH_MANIFEST.json),
[CRC/full-object-SHA receipt](delivery/EVIDENCE_ARCHIVE_RECEIPT.json).
Content-addressed objects retain every original included path, including identical
bases; no original raw payload is omitted. The manifest and packaging source
explain excluded circular archive/receipt metadata. Binary/source freezes,
all55 raw attempts, every audit, validation sidecar and report/graph are included.

Post-campaign release library4089passed/0failed/117ignored. All48site tests pass;
generated leaderboard/browser/panel index are current. New native coverage,
recursive A/A and raw-SHA tamper tests pass. Original source/binary/v1 seals
pass unchanged. Targeted unchanged curve-name lint passes all changed/new
Markdown/HTML; earlier unrelated sparse full-lint access failures are recorded,
with no lint/CI guard altered. [Source/PDF review](checks/final/PDF_AND_SOURCE_REVIEW.md),
[canonical diagnostic review](checks/final/CANONICAL_DIAGNOSTICS.md).

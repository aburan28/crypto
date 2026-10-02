# Fused sigma schedule: measured result

Decision: **promote as compatible engineering.** The fused schedule removes
one coordinate load pass from every complete sigma update and preserves the
walk, checkpoint-visible point state and corpus. It improves the counter-free
headline benchmark from 15.015004 to **15.436677 billion complete scalar
updates/s** on one RTX PRO 6000. The active 26 B/s goal remains unachieved.

## Single table

Every timing sample completes 403,726,925,824 updates. Medians are over five
alternating A/B pairs after four excluded warmups. A/A is the same control
binary run in both positions.

| workload | control B/s | fused B/s | paired median | paired minimum | A/A max drift | replay/corpus | class |
|---|---:|---:|---:|---:|---:|---|---|
| counter-free headline (`WITNESS=0`) | 15.015004 | **15.436677** | **1.028257** | 1.027890 | 0.0861% | 300/300 each, zero drops, identical 1,709,477-record v1 corpus | engineering; selected |
| witness-carrying (`WITNESS=1`) | 9.347867 | **11.800903** | **1.262291** | 1.261768 | 0.0574% | 300/300 each, zero drops, identical 1,709,477-record v2 corpus | engineering; selected for this mode |
| active objective | | **26.000000** | headline / goal = 0.593718 | | | not met | boundary |

The five headline ratios are 1.028140, 1.028257, 1.028864, 1.028329 and
1.027890. The five witness ratios are 1.262291, 1.262438, 1.262145, 1.262538
and 1.261768. Every pair favours fusion, and both A/A panels are well inside
the frozen 1% noise gate.

The counter-free row is the headline because it matches the existing campaign
and throughput objective. It is 3.917% above the prior 14.8547595 B/s verified
sigma result. Reaching 26 B/s from the new best still requires 1.6843x.

## Refreshed inline control

The fused comparison deliberately holds `PACKED_INLINE_POLY=3` in both arms.
A separate current-source five-pair panel re-establishes that prerequisite
against out-of-line polynomial helpers: **14.941312 versus 14.534256 B/s**,
with paired ratios 1.025384/1.028820/1.028007/1.028821/1.029012 and a paired
geometric mean of 1.028008. Its 95% paired log-ratio interval is
`[1.026123, 1.029896]`. Both arms replayed 300/300 with zero drops, passed
arithmetic/storage/shared-sigma and bidirectional checkpoint gates, and
produced identical 29,779-record v1 corpora.

That producer used the repository's legacy format-aware comparison script.
Before incorporating its result here, the native C++ comparator in this
directory independently reopened both raw corpora, confirmed headerless-v1
framing and byte-identical sorted records, and reproduced canonical SHA-256
`17e29b9695e175f306b0c76136b7331601596c98b4a9e1c27311f03c6b4a635d`.
The successful archive and the preserved first-attempt producer failure are
bound in `results/inline-refresh-artifact.json`.

## Closed geometry and direct-map follow-ups

The missing matched B16 block-width comparison is now closed.  On one current-
main allocation, B16/T512/min1 beat B16/T256/min2 in all five pairs, but only
by a 1.006030 median ratio (minimum 1.005781) against 0.0833% maximum A/A
drift.  That misses the frozen 1.01 selection threshold, so the native fused
preset remains **B16/T256/min2**.  Both arms replayed 300/300 with zero drops,
produced identical sorted 1,711,916-record corpora, and passed cross-geometry
checkpoint replay.  Protocol, raw manifests and the corrected native audit are
in [`b16-geometry/`](b16-geometry/).

The larger direct polynomial-basis `I + sigma^j` experiment is also terminal.
Its exact three-bit shared table cut fixed-`j` host instructions by 30.2%, but
on the RTX PRO 6000 it reduced matched throughput from 15.893066 to 13.143106
B/s: median ratio 0.826916, with every pair negative and 0.0371% A/A drift.
All map, replay, corpus and checkpoint gates passed.  The default remains the
composed shared-sigma network; terminal evidence is in
[`direct-map-rejection/`](direct-map-rejection/).

## Correctness and equivalence

Both protocols use seven odd 95-step launches at DP weight 48 before timing.
Each arm replays 300 reports against the host reference with zero drops. The
counter-free arms contain 1,709,477 headerless 32-byte v1 records and share the
sorted-payload SHA-256
`7e0e8c90a7a1256c4f2c7f5ad4f23d5334de025a82590639f25516c0be25de92`.
The witness arms contain 1,709,477 framed 72-byte v2 records and share
`9da10fbdd66dc175482e4e100fea68f541020881fbdea9bbc1bf7a4c9087e5c7`.

The independent witness audit parsed every record: all 1,709,477 satisfy
`sum(branch_counts) == iters`; the 5,300 zero-count records are legitimate
zero-step distinguished starts. The other records account for 395,204,474
reported trail steps, with maximum trail length 664.

An independent native model uses the repository's actual packed GF(2^131)
arithmetic and compares the ordinary schedule with the alternating fused
schedule for 1, 2, 3, 4, 7, 16 and 95 steps over seven launch boundaries. It
passes for both global and shared sigma masks, including identical per-slot
DP, guard and witness observation times. Inactive partial-block threads still
participate in the shared-mask initialization barrier before returning.

## Resources and scope

The headline control uses 91 registers/thread; the fused candidate uses 126.
Both report zero local bytes and 1,792 shared bytes per block and retain two
256-thread blocks per SM. The witness candidate uses 120 registers.

This is one-GPU kernel engineering in the existing sigma walk. It changes no
collision count, automorphism quotient, recovery formula or ECDLP exponent. No
search, solver or collision recovery was run.

Independent audits:

- `results/headline-independent-audit.json`, SHA-256
  `1b1a46ca52aa14ef86b96fc9b765c9c10e2815c1fe9bcad867422785ab7976d1`
- `results/witness-independent-audit.json`, SHA-256
  `993463aeec003de0475fd2dc552ffa6757715abd7b1992dc6022930bc61ed076`

Raw archives remain retrievable through the Modal volume tokens in
`results/headline-artifact.json` and `results/witness-artifact.json`; their
SHA-256 values are `40881036…97e1` and `abb45744…69ac`. The archives record
binary hashes but do not embed executable bytes. Their in-job manifests hashed
`job.log` before the wrapper finished populating it and omit the later
`exit-code` file; the committed post-run file manifests correct that packaging
race. The headline artifact's preflight text says “v2,” a label typo: the
protocol, native comparator and independent parse all establish v1 framing.

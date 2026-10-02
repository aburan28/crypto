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

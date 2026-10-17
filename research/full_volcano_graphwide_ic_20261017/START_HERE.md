# Volcano experiment downloads, 17 October 2026

[Open the Drive folder](https://drive.google.com/drive/folders/109DpcQpX7ATz6GeI8iSWAExGMWBoR-H_).
These are frozen snapshots. Complete graph enumeration and completed ECDLP
measurements have different coverage.

| Artifact | What it contains |
| --- | --- |
| [Latest complete-graph explorer](https://drive.google.com/file/d/1JT72Pxe6cwJwlr7WBfvbu5eC4Q1fNp7E/view) | All 24,823 vertices / 69,300 selected-degree maps; frozen 2,935 audited-vertex overlay. Download HTML and open locally. Not full measurement coverage. |
| [Original C prefix dataset](https://drive.google.com/file/d/1Xi7mjLCskYEauXlwzAfQQtOoHw6waop-/view) | 1,869 known C vertices; 3,738 verified finite-field relation targets, exact re-execution, binaries, native statistics, local certificates, calibration graphs and report. Largest prime subgroup 7. NOT full C class. |
| [Original C report](https://drive.google.com/file/d/1jqNQUqu2uxOU3n7NOn1jOYFQaAD9We4L/view) | Four inspected pages. Local exact conductors: 517 vertices / 24 strata, 2,985,984–967,458,816; 1,352 unknown. Source conductor 967,458,816. |
| [C per-target operations CSV](https://drive.google.com/file/d/1B4XZNFrlzgIzZT1vI7gq6kJs1MfWFFQf/view) | All 3,738 observations with phase/native counters, unqualified elapsed intervals and local conductor annotations. |
| [Repaired runner QA](https://drive.google.com/file/d/12Q0z7VUCWJuCgY5ofbW5Dof3V15FZ4Th/view) | New control and pilot: 60 attempts (40 verified, 18 exhausted, 2 timeouts); 58 completed baseline matches, two censored comparisons. Whole-queue execution is separate. |
| [Stopped release checkpoint](https://drive.google.com/file/d/1vVW3fCYqdnl3Ue2ef1FckbSTODbg-MTs/view) | 70 audited vertices / all 51 strata, 1,120 accepted attempts; extra failed 16-record capsule preserved and excluded from accepted means. Do not resume its old runner. |
| [Historical dev99 graph explorer](https://drive.google.com/file/d/13AQHFuLnePIfbbS9XbaXNueJY-JMyKCE/view) | All 24,823 vertices, 51 conductor strata and 69,300 recorded directed kernel maps. Download the HTML and open it; Drive preview does not run it. |
| [99-vertex measured archive](https://drive.google.com/file/d/1hLVWVad2SrgHaBn6U5iSyqnKbx5cAkeY/view) | Complete three graph databases, exact frozen dev binaries/source, full 24,823-row queue, execution ledger, 99 sealed sessions, replay receipts, native statistics and report. |
| [Three-page report](https://drive.google.com/file/d/1h_G5uJZYYqVzBL3_kynRdYZvE8vjJu9z/view) | Measured coverage and yield/cost summaries over all 51 strata; every page inspected. |
| [Detailed statistics CSV](https://drive.google.com/file/d/1OWEkYFgWuK_wgMu0auEHsYxJ2VDugFIA/view) | Per-graph/stratum/arm/outcome min, max, arithmetic mean and sample count, including phase/native-operation counters. |
| [Native pricing band](https://drive.google.com/file/d/1NV-4hpuoG23Xf0TitT8yUBPGAnRJCoj2/view) | All 792 curve/arm rows at common unpriced-operation weights 0/0.1/1 GAE; incomplete/unaccounted cases remain explicit. |
| [Separate optimized control/pilot](https://drive.google.com/file/d/1lJt6NhNwXVpvdSBoXqCXhoi-flqEzwmF/view) | Actual release runs, binaries, exact dev baseline capsules and native comparison receipts. This is not the optimized whole-queue survey. |

The three source models are `icv1-fp57-t536870914-6cabb0c9`,
`icv1-fp57-tm536870941-e6d0c766` and `icv1-fp59-t1073741954-4dd0b1ef`.
Their maximum endomorphism conductors are 8192, 6561 and 1728, respectively.
Their largest prime subgroups have orders 61751, 19768889568671 and 4253;
field size is not subgroup security size. Every vertex uses its certified
actual largest-prime subgroup, not a substituted smaller subgroup.

The stopped dev checkpoint has **99 measured vertices across all 51 strata**:
1,584 attempts, 1,065 verified, 422 exhausted and 97 timeouts. It compares
projected factor bases 16/64, MITM/subtraction, a raw-base control and matched
rho/BSGS, using two public targets, one round and bounded trials/watchdogs.
Every completed successful/exhausted deterministic report has a replay receipt.
These are same-implementation replays, not independent validation.

It found 57,176 relations from 2,223,204 source trials and 2,075,145 actual
decomposition calls. Repeated queries skipped before the oracle explain the
different denominators; summary v2 exposes both. Timeout counters are unknown,
not zero. The old first-12 report's "decompositions" label meant source trials;
historical records/reports/uploads remain unchanged. Do not pool different
subgroups, build profiles or repeated control/pilot targets as one population.

The optimized control/pilot have 60 attempts: 40 verified, 18 exhausted and
two timeouts. All 56 comparable completed deterministic profile pairs match
exactly; four baseline-censored pairs have no final counts to compare.
The first optimized survey halted on a failed audit after 70 accepted vertices;
its late-message watchdog recording defect is repaired without relabeling old
records or weakening the auditor. A separately frozen repaired whole queue
uses a new binary/ledger and 24-hour/2-GiB guard; actual launch/progress require
its separate launch receipt and local ledger. Its results are separate from this dev snapshot. No runtime
speedup, fastest-method or asymptotic claim is established by these runs.

The graph CLI is `crypto isogeny class`. For a supported ordinary prime-field
curve JSON, use `plan --curve curve.json`, then
`build --curve curve.json --out graph.sqlite --steps all`. Inspect with
`vertex --input graph.sqlite --id N` and `status --input graph.sqlite`.
The native [`CLI documentation`](https://github.com/aburan28/crypto/blob/codex/full-isogeny-class-cli-20261016/tools/volcano-cli/README.md)
defines backend limits, resume, map verification and the separate one-target
BSGS command. Graph generation does not itself run graph-wide index calculus.

Graph completeness is for Fp-isomorphism vertices and the declared finite
degree sets, not other primes or infinitely many compositions. The original
A–E graphs/campaign, full eight-target/five-round/million-trial matrix, exact
per-vertex representation registration, transported targets, composite PH and
independent replay remain open. The C k=1 diagnostic above is implemented and measured; finite-field IC for
the other actual embedding degrees remains unimplemented. The k=2 MOV teaching
example is not those backends. The C method is projected-prime subgroup
relations, not NFS or integer-smoothness IC.

Archive extraction, manifest checks and repacking validations are scoped
individually in the delivery receipts. The new 2,935-vertex and five-prefix
archives pass extraction and every manifest entry; their outer repacking was
not rechecked. Drive metadata confirms names, sizes and parent folder; remote
whole-file hashes were not checked. Delivery JSON files retain exact byte
counts and SHA-256s.


Original A–E total 7,833,658,967 expected vertices and 515 strata. Historical
B is separately complete; A/C/D/E remain incomplete. Google Drive deliveries
work. This runtime's S3 probe found no AWS identity, so no S3 upload is claimed.
The repository's reports-to-S3 workflow runs after merge to main; it does not
provide a larger random-access SQLite working volume. A connected existing
worker/volume and S3 identity remain required for those large runs.


Historical launch: repaired local worker **PID 79776 launched 2026-10-17 06:58:03 UTC**, with
24-hour / 2-GiB guards. At 06:59:13 UTC it has audited vertices in all three
classes and an active job. This is an actual launch, not full-queue completion.
That process was absent on recovery; this is not its current running status.
Later progress is in its local ledger and `execution/repaired-continuation/`.
It will pause on guards/failures; no unattended upload is assumed. The frozen
Drive snapshots do not track live progress.


Latest checkpoint, 17 October 2026 14:50 UTC: the earlier process was absent
after the environment returned; its PID receipt is historical. A preserved
interrupted attempt and **2,935 audited vertices / 46,960 records** are in the
[new checkpoint archive](https://drive.google.com/file/d/1IoPt5C-EYlo488pPVWPx-qpDHo7t1GwY/view).
[Latest all-arm min/max/mean CSV](https://drive.google.com/file/d/1AK_WDN-qMx_2FGkDkgfIqH5PyP6gjUus/view).
There are 30,270 verified recoveries, 14,626 exhausted budgets and 2,064
timeouts, with samples from every one of the 51 strata; 1,530,557 relations.
The [report](https://drive.google.com/file/d/1fiyasnknjZ1T0lE8z6yHi7AIhMdUysnA/view)
is inspected. New runtime epochs and resumes have separate receipts.

[Original A–E local conductors and prefixes](https://drive.google.com/file/d/1zeGSq3GWne9yTH0R37EZ3YI89xG_yR2L/view)
now include source conductors 268,435,456 / 129,140,163 / 967,458,816 /
544,195,584 / 644,972,544 and per-vertex bounds, with unknowns retained.
[One-page report](https://drive.google.com/file/d/17jC_A1uSRHtJc5c8Sy9AkH9jt6UJTbBX/view).
This is local certification, not complete original class populations. S3
access was rechecked and still lacks a runtime AWS identity.

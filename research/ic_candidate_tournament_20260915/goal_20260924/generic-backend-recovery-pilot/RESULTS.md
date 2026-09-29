# Disclosed recovery pilot: F5 complete, combined F4/F5-and-SAT gate failed

The [preregistered panel](panel.json) ran all twenty jobs once on five
previously disclosed public points. One n17a1 F5 worker recovered and
independently replayed the logarithm `40605` through full-rank relation
collection, final linear algebra and target descent. Neither n17a1 SAT arm
returned a report before its 300-second cap. Therefore the registered
**combined admission condition was not met**. This is a local scientific
admission result, not a fresh-target, incumbent, rho or speedup comparison.

The [complete raw bundle](evidence.tar.gz) retains the exact jobs, worker
stdout/stderr, process records, audit receipts, frozen panel and runner,
source/dependency manifest, source snapshot, build policy, executable and
progress/summary records. Its SHA-256 is in [EVIDENCE.sha256](EVIDENCE.sha256).
`test_generic_backend_recovery_evidence.py` verifies the archive, source
snapshot and executable hashes, replays the independent admission of every
worker report, checks all twenty fixed schedule rows and confirms paired
ordinary-query sequences. Timed-out workers have no report or inferred PDP
outcome. The measured build uses pinned source commit
`765c3c5f19032bd852163805f257c56babef2040`, source-manifest hash
`c64e4b3102bface63a2305efbff4bd85810cc112cb43546da2992ba48e9e85b7`
and worker hash
`94caf3d67e57dde09488763ec19e791aca54bbb35f1968c73e42d37a210b436a`.

The dimension-six standard-subspace bases were independently enumerated and
subgroup-filtered before the run. These counts are the controlled fixture
inventory; a timed-out arm did **not** independently report its own base.

| Curve cell | Actual usable points before folding | Sign/Frobenius columns |
| --- | ---: | ---: |
| n17a1 | 62 | 29 |
| n19a0 | 62 | 27 |
| n23a0 | 72 | 33 |
| n23a1 | 52 | 23 |
| n31a0 | 66 | 27 |

`W/I` counts verified ordinary-query witnesses and budget-incomplete PDP
attempts. Every reported attempt, including an incomplete one, is retained.
`—` means no auditable report, not zero yield. Process wall is a local
diagnostic that includes target-independent preparation; online wall is
reported only for the complete, audited one-target solve.

| Cell | Solver | Process result | W/I | Final rank/columns | Process wall (s) | Verified online wall (s) |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| n17a1 | F4 | timeout, no report | — | —/29 | 300.069 | — |
| n17a1 | F5 | **verified complete** | 29/75 | **29/29** | 273.978 | **15.236** |
| n17a1 | SAT native-XOR | timeout, no report | — | —/29 | 300.040 | — |
| n17a1 | SAT CNF | timeout, no report | — | —/29 | 300.019 | — |
| n19a0 | F4 | timeout, no report | — | —/27 | 300.080 | — |
| n19a0 | F5 | bounded incomplete | 2/22 | 2/27 | 206.144 | — |
| n19a0 | SAT native-XOR | bounded incomplete | 0/24 | 0/27 | 112.396 | — |
| n19a0 | SAT CNF | bounded incomplete | 0/24 | 0/27 | 184.625 | — |
| n23a0 | F4 | bounded incomplete | 0/4 | 0/33 | 62.582 | — |
| n23a0 | F5 | bounded incomplete | 0/4 | 0/33 | 31.758 | — |
| n23a0 | SAT native-XOR | bounded incomplete | 0/4 | 0/33 | 11.442 | — |
| n23a0 | SAT CNF | bounded incomplete | 0/4 | 0/33 | 16.771 | — |
| n23a1 | F4 | bounded incomplete | 0/4 | 0/23 | 61.813 | — |
| n23a1 | F5 | bounded incomplete | 0/4 | 0/23 | 31.515 | — |
| n23a1 | SAT native-XOR | bounded incomplete | 0/4 | 0/23 | 10.493 | — |
| n23a1 | SAT CNF | bounded incomplete | 0/4 | 0/23 | 24.844 | — |
| n31a0 | F4 | bounded incomplete | 0/1 | 0/27 | 70.799 | — |
| n31a0 | F5 | bounded incomplete | 0/1 | 0/27 | 27.686 | — |
| n31a0 | SAT native-XOR | bounded incomplete | 0/1 | 0/27 | 2.651 | — |
| n31a0 | SAT CNF | bounded incomplete | 0/1 | 0/27 | 5.307 | — |

All sixteen exited reports passed the source-bound build, actual base,
ordinary-query/group-relation, observed solver-dispatch, matrix/rank and
exclusive-phase audits. The n19–n31 arms with reports used the same ordered
ordinary queries within each cell. Four processes timed out without stdout;
their attempted query count and natural yield are unknown. There was no OOM
or unsupported-solver report. No timed-out job was retried.

The n17a1 F5 run used 104 ordinary collection attempts to obtain 29 verified,
linearly novel relations, with 75 bounded-incomplete attempts. Its retained
certificate records 62 usable base points, 29 folded columns, rank 29,
29 verified relations, one certified descent and one scalar replay. The
exclusive worker ledger closes at 273.978 s: 258.655 s of collection PDP
(including failed attempts), 15.236 s of target descent work, and the
remaining setup, query, base, matrix and checks. The online interval is
15.236145209 s, of which 15.236103625 s is target PDP; it starts after
reusable factor-base logs are ready and includes the worker's scalar replay.
The local process wall, memory sample and online interval have one observation
and no calibrated Linux instruction count. All incomplete runs retain their
observed phase costs but have missing final LA/descent costs, so their full
cost and speedup remain unknown.

For an exploratory description of observed witness frequency, a nominal 95%
Wilson interval gives 29/104 = 27.9% (20.2–37.2%) for n17a1 F5,
2/24 = 8.3% (2.3–25.8%) for n19a0 F5 and 0/24 (upper 13.8%) for each
n19a0 SAT arm. The zero-yield n23 cells have only four trials per arm (upper
49.0%); n31 has one (upper 79.3%). These intervals assume independent
exchangeable queries, while the stream was fixed and n17 F5 stopped at full
rank. They describe this disclosed-point pilot only; they do not measure
between-target variation or establish fresh-target yield. The n19 paired
process times have no repetition-based uncertainty estimate. A quicker
zero-yield SAT process is not a faster IC solve.

No comparison to incumbent or paired rho was run. The one-target online
speedup, operation-normalized `S`, global tournament rank and promotion status
are unknown. The next admission step is to investigate why the present SAT
encodings exhaust their bounds, preregister a changed SAT backend or encoding
with a full source/build manifest, and obtain an independently verified
complete disclosed-point SAT solve before spending fresh paired targets.
The separate live generic-v2 seed `2026092902` remains a distinct frozen
campaign and must be audited at completion without redispatch.

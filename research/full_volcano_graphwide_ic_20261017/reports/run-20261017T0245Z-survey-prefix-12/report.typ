#set page(paper: "a4", margin: 13mm)
#set text(size: 8.5pt)
= All-vertex survey: first 12 audited strata
Run 20261017-survey-prefix-12. Frozen checkpoint, not completion.

#image("strata-yield.svg", width: 100%)

*Native recorded coverage:* 12 audited vertices in 12 distinct graph/stratum
pairs. Four vertices in each graph. The immutable queue contains all 24,823
vertices and 51 strata; the first phase visits one vertex per stratum before
all remaining vertices. Pilot/root reruns use the same targets and seeds;
they are repeated observations, not additional independent targets.

*Observed outcomes:* 192 measured attempts: 119 verified, 58 exhausted,
15 watchdog timeouts. All 12 standard audits pass. Each successful deterministic
record replays, and each exhausted deterministic record passes the additional
native replay receipt. No wrong answers were recorded. Timeouts have unknown
final relation and operation counters. Their elapsed observations are retained.

*Recorded direct-IC counters:* 301,715 attempted decompositions, 6,368 relations
found across the 120 direct-IC attempts with final reports. Numerator and
denominator include exhausted reports and the raw-base control. They are not
120 successful DLPs or independent identically distributed trials. Phase,
rank, dependency, failure, memory and unpriced-work records remain available.

*Subset below:* size-64 projected MITM only. Cold GAE are lower bounds with
unpriced native work; no fully priced cost or runtime ranking follows.
Each row is two measured one-target attempts on one exact vertex.

#set text(size: 7.5pt)
#table(columns: (.8fr,.5fr,1fr,.3fr,1fr,1fr,1fr,1fr), inset: 3pt,
  [Graph key], [f], [Status], [n], [GAE min], [GAE mean], [GAE max], [Mean yield],
[two], [1], [verified], [2], [19,530], [19,970.0], [20,410], [0.119769],
[two], [2], [verified], [2], [19,505], [19,546.5], [19,588], [0.114059],
[two], [4], [verified], [2], [19,466], [19,531.0], [19,596], [0.123116],
[two], [8], [verified], [2], [19,467], [19,537.5], [19,608], [0.128813],
[three], [1], [exhausted], [2], [20,239], [20,252.0], [20,265], [0.000000],
[three], [3], [exhausted], [2], [20,245], [20,268.0], [20,291], [0.000000],
[three], [9], [exhausted], [2], [20,249], [20,251.0], [20,253], [0.000000],
[three], [27], [exhausted], [2], [20,249], [20,259.5], [20,270], [0.000000],
[mixed], [1], [verified], [2], [20,463], [20,467.5], [20,472], [0.840387],
[mixed], [2], [verified], [2], [20,216], [20,228.5], [20,241], [0.821442],
[mixed], [3], [verified], [2], [20,281], [20,438.0], [20,595], [0.856216],
[mixed], [4], [verified], [2], [20,095], [20,123.0], [20,151], [0.840969],)
#set text(size: 8.5pt)

*Conditions:* actual largest-prime subgroups; projected sizes 16/64 with
MITM/subtraction, unchanged raw-base MITM, rho.negation, textbook BSGS and rho
A/A. Two public targets, one round, no warmup, 4,096 trials, 3-second watchdog.
Prime-subgroup membership and all recovered scalars are checked. Every run
charges its cold solver setup; no reusable logs are assumed free. Ordinary
finite-field IC remains unsupported at the actual recorded embedding degrees.

*Host:* Intel Xeon Platinum 8573C, Linux x86-64, four-CPU quota, 16 GiB RAM.
Unoptimized ecbench SHA-256 1869ca2e7be271794df8128c08a6dca0201038e2c81a3e0ce6434f7ec0d8379d;
native driver source commit ee54ad603. All records earn L0; raw solve/process
times are descriptive, not admitted performance. Short driver rejection checks
ran during this counts-only block; no build or test suite ran during it.
The observed invocation took 272 elapsed seconds, including runner overhead
and audits. This is a local scheduling observation, not an all-vertex forecast.

*Evidence:* execution/campaigns/survey-prefix-12, immutable full queue in the
delivery archive, campaign.json hashes, SQLite checkpoint, standard audits
and exhaustion-replay receipts. The native summary v2 includes queued/audited
coverage by conductor and explicit sample counts for every min/max/mean.
Common 0/0.1/1 GAE pricing-band output is separate. Same-implementation replay
is not independent implementation validation. A second writer and changed
driver were rejected; zero-work resume preserved the 12 audited sessions.

The original full matrix (million trials, eight targets, five rounds, warmup,
four base sizes), original A--E full graphs, exact representation registration,
transported targets, full composite PH, qualified runtimes, finite-field IC
and independent implementation replay remain open. A real queue is not a
completed solve. Canonical best-method/frontier claims are unchanged.

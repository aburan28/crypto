# SAT/PDP matched-pipeline workstream

The first preregistered protocol is [`PROTOCOL.md`](PROTOCOL.md). Its
[`pilot`](sessions/pilot) closed with 28 records: 24 verified and four
`ic-buchberger` timeouts at the fixed 60-second per-run limit. The
strong-rho, meet-in-the-middle, Frobenius meet-in-the-middle, and both
SAT encodings returned verified answers on both public targets in
both rounds. The timeout stops the original eight-target main panel
under its stated rule; [`specs/main.json`](specs/main.json) remains
frozen and unrun.

The first pilot also exposed a cost-replay defect in the historical
`ic.pipeline` method. Its
[`failed audit`](pilot.audit.failed.json) checked every seal and identity,
then replayed 11 of 12 measured verified runs identically. One SAT
native-XOR replay had the same answer and integer counters but differed
in the low bits of `total_gae` and phase GAE. Receipt SHA-256:
`7a25f644ddc2fc8137d2a4e2c1fbd54993798eeb3b543c2d28a2effa3b6feef1`.
That failed receipt is retained; the original pilot is not used for a
speed ratio.

[`PROTOCOL-2.md`](PROTOCOL-2.md) is a new, narrower preregistration.
It drops the timed-out Buchberger arm and uses the versioned
`ic.pipeline_counted` method, which reconstructs charged relation
cost from counts instead of subtracting a wall-derived floating charge.
The method keeps SAT conflicts and other units without a pinned
conversion explicitly unpriced. Its pilot and holdout panel have
separate frozen specs. The holdout panel's result is pending; the
first pilot cannot stand in for it.

The counted-method [`pilot`](sessions/counted-pilot) completed with 24
of 24 verified executions. Its [`full-replay receipt`](counted-pilot.audit.json)
passed with all 12 measured records identical and zero problems;
receipt SHA-256
`090ae757b4899c236a01a1515ec445f802e097e94ee9be4590008f8aa3bfca21`.
It exposed a missed `PROTOCOL-2.md` stop gate: the two measured CNF
runs had 56 per-call conflict-budget exits (110 across warm-ups and
measurements). Verified final answers do not erase those exits. The
six-arm main panel was started in error and stopped; its 208-record
stale attempt and 80-record cleanly interrupted restart are retained
with hashes and an interrupted-session replay receipt in
[`INTERRUPTED.md`](INTERRUPTED.md). Neither enters a comparison.

[`PROTOCOL-3.md`](PROTOCOL-3.md) was frozen before the follow-on run.
It excludes only the budget-hitting CNF arm and keeps the remaining
methods, public targets, factor base, seeds and budgets fixed. Its
[`count-qualified-main`](sessions/count-qualified-main) session completed
with 240/240 verified executions, 40 measured runs for each of five
arms, zero trial or solver budget exits and zero timeouts. The
[`full-replay audit`](count-qualified-main.audit.json) reports 200/200
measured runs identical and zero problems. Receipt SHA-256:
`a03b673ed8f2ee75bdf53c4c8e13c385b0f2c506361d8d28158f9b579324b1e4`.

| arm | mean cold S | S / matched strong rho | S / generic floor | verified | class |
|---|---:|---:|---:|---:|---|
| strong signed-Frobenius rho | 3.998 | 1.000 | 18.598 | 40/40 | reference; native work unpriced |
| rho A/A control | 3.998 | 1.000 | 18.598 | 40/40 | control; counted work identical |
| IC MITM | 387.277 | 96.880 | 1801.777 | 40/40 | baseline |
| IC Frobenius MITM | 17.073 | 4.271 | 79.430 | 40/40 | engineering |
| IC native-XOR SAT | 3.413 **lower bound** | 0.854 **diagnostic** | 15.879 **lower bound** | 40/40 | accounting incomplete |

The <a href="count-qualified-main.table.md">native ecbench table</a> and
the seven saved <a href="sessions/count-qualified-main/comparisons/">paired comparisons</a>
carry the full intervals and per-workload rows. All operation ratios to
rho are diagnostics because its native table and canonicalisation work
is unpriced; the SAT row additionally omits solver work. Dividing
lower bounds gives no bound on the true cost ratio. In particular,
the SAT row's counted `S` below rho is **not** a full-pipeline or
runtime speedup.

The <a href="count-qualified-main.resources.json">native resource vector</a>
records 2,755,966 SAT conflicts in 749 solver calls, zero conflict-budget
exits, and 8,560 KiB peak child RSS for native-XOR SAT. All 200 measured
runs earned L0 on macOS; none met the required L2 wall gate, and this
host exposed no solve PMU counts. The A/A operation ratio is exactly
1.000. Its L0 wall median is descriptive only. Every session record
contains its cold phase costs and one-target online window. The result
is an accounting demonstration on one small Koblitz curve, not a
cross-size claim or an m = 83 confidence result.

The same resource report separates 40-run cold charged phase totals:
strong rho used 30,722.155 GAE in set-up and 10,228.202 in search;
IC MITM used 3,933,840 in oracle set-up and 30,382.734 in relations;
Frobenius MITM used 141,480 and 30,383.709; SAT used zero in oracle
set-up and 31,929.963 in charged relations **plus** its unpriced
conflicts. The exclusive online-window sums are 1,222,546 ns for rho,
8,066,083 ns for MITM, 6,858,326 ns for Frobenius MITM and
41,620,896,790 ns for SAT. These L0 times are phase diagnostics only.
At this small subgroup size the strong rho's charged set-up dominates
its own `S`, so the panel cannot support a scaling inference.

The frozen session names source commit
`b83ff0902ee0d7d1d7fbed07f92ce441ee974f92` (clean), binary SHA-256
`4c6f3b6d45f988ad134c6055600014e236858fa674b03ce7785a8b320b325869`,
spec SHA-256
`ab496597915646afba9b7905a2046ea896c6317cd0c69737b3193773e0751b00`,
records SHA-256
`8f499a1d6f98164425bdc90d1903d935ca36686f2a36caacf2180fcee45d2f59`,
and host class `ECBENV2hd82681268e96`. A second architecture's CI
replay, real two-socket NUMA/PMU measurements, calibrated native-work
prices, and the wide m = 83 IC/strong-rho pair remain pending gates.

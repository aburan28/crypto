#set page(paper:"a4",margin:13mm)
#set text(size:9pt)
= Watchdog repair: control and pilot evidence
17 October 2026. New campaign data, separate from historical dev/release runs.
#image("qa-outcomes.svg",width:100%)

The earlier release survey halted after 70 audited vertices on a failed
standard audit. One late child message carried a checked scalar while the
parent labeled the attempt timeout. Its failed 16-record capsule is retained;
no historical record or auditor rule was changed.

The repaired runner makes watchdog censoring authoritative. Timed-out attempts
have no accepted scalar, success check, final cost/counters or solve timing.
Any parsed late message remains in the sealed diagnostic error. Solver
arithmetic, public targets, factor bases and budgets are unchanged.

Fresh small control: 12 attempts. Fresh three-root pilot: 48 attempts. Every
session passes the standard audit; successful and exhausted deterministic
reports reproduce. Native comparison to the earlier release baseline checks
method identities, exact outcomes, scalars, GAE, phases, field operations,
counters and factor bases. It finds 12 matching completed control pairs and
46 matching completed pilot pairs. Two censored pilot pairs have no completed
reports to compare. Elapsed values are excluded. Same-implementation replay
and profile matching are not independent implementation validation.

The frozen runner hash is
#text(size:7pt)[35d89073f4d666b094c213a1d31b0c819b63600123e067c111af3c947d7609a9].
The driver hash is
#text(size:7pt)[b8b489c7225c64ed492ffce39286c995f5ca820960956cf706c0d3110f5bdd2d].
Clean source 2bf7700b7 contains the repair. The native build source/protocol
and input/binary hashes are frozen in each campaign. The separate full serial
release library check passes 4,236 tests, zero failures, 117 existing ignored.
The previously failing parallel Joux--Vitse case passes serially; its parallel
cause remains unresolved. No solver budget or exact-answer assertion changed.

Host: Xeon Platinum 8573C, Linux x86-64, four-CPU quota, 16 GiB memory limit.
Release optimization; wall observations remain unqualified. No external
comparison, fastest-method claim or speedup is made. The methods retain actual
largest prime orders 61751, 19768889568671 and 4253 in the pilot; the small
control has order 7 over p=113. No smaller subgroup replaced the 45-bit case.

Passing these gates permits a separately frozen 24,823-vertex bounded survey,
not the full eight-target/five-round/million-trial matrix. Planning is not
execution. Its actual launch and current coverage require the separate launch
receipt and live ledger. Original A--E completion and large-extension
finite-field IC remain open. Source protocols, raw capsules and native
comparison receipts are in the accompanying archive and study execution tree.

#pagebreak()
= Native per-arm GAE distributions
#set text(size:6.7pt)
#table(columns:(.65fr,1.5fr,.9fr,.35fr,.35fr,.85fr,.9fr,.85fr),inset:2pt,
 [Graph], [Arm], [Outcome], [Attempts], [Cost n], [Min GAE], [Mean GAE], [Max GAE],
[control], [bsgs], [verified], [2], [2], [4], [4.00], [4],
[control], [projected-2-mitm], [verified], [2], [2], [120], [123.50], [127],
[control], [projected-2-subtract], [verified], [2], [2], [111], [115.00], [119],
[control], [raw-control], [exhausted], [2], [2], [1,327], [1,335.50], [1,344],
[control], [rho], [exhausted], [1], [1], [8,013], [8,013.00], [8,013],
[control], [rho], [verified], [1], [1], [23], [23.00], [23],
[control], [rho-aa], [exhausted], [1], [1], [8,013], [8,013.00], [8,013],
[control], [rho-aa], [verified], [1], [1], [23], [23.00], [23],
[mixed], [bsgs], [verified], [2], [2], [75], [104.00], [133],
[mixed], [projected-16-mitm], [verified], [2], [2], [4,242], [4,419.00], [4,596],
[mixed], [projected-16-subtract], [verified], [2], [2], [10,148], [12,328.00], [14,508],
[mixed], [projected-64-mitm], [verified], [2], [2], [20,463], [20,467.50], [20,472],
[mixed], [projected-64-subtract], [verified], [2], [2], [22,087], [22,588.00], [23,089],
[mixed], [raw-control], [exhausted], [2], [2], [13,361], [13,383.00], [13,405],
[mixed], [rho], [verified], [2], [2], [138], [189.50], [241],
[mixed], [rho-aa], [verified], [2], [2], [138], [189.50], [241],
[three], [bsgs], [timeout], [2], [0], [missing], [missing], [missing],
[three], [projected-16-mitm], [exhausted], [2], [2], [8,191], [8,204.00], [8,217],
[three], [projected-16-subtract], [exhausted], [2], [2], [138,735], [138,748.00], [138,761],
[three], [projected-64-mitm], [exhausted], [2], [2], [20,239], [20,252.00], [20,265],
[three], [projected-64-subtract], [exhausted], [2], [2], [536,271], [536,284.00], [536,297],
[three], [raw-control], [exhausted], [2], [2], [14,479], [14,492.00], [14,505],
[three], [rho], [verified], [2], [2], [2,005,511], [3,519,227.50], [5,032,944],
[three], [rho-aa], [verified], [2], [2], [2,005,511], [3,519,227.50], [5,032,944],
[two], [bsgs], [verified], [2], [2], [266], [276.50], [287],
[two], [projected-16-mitm], [verified], [2], [2], [6,493], [6,927.00], [7,361],
[two], [projected-16-subtract], [verified], [2], [2], [76,008], [90,295.50], [104,583],
[two], [projected-64-mitm], [verified], [2], [2], [19,530], [19,970.00], [20,410],
[two], [projected-64-subtract], [verified], [2], [2], [130,718], [160,668.50], [190,619],
[two], [raw-control], [exhausted], [2], [2], [13,726], [13,750.00], [13,774],
[two], [rho], [verified], [2], [2], [355], [360.50], [366],
[two], [rho-aa], [verified], [2], [2], [355], [360.50], [366],
)
#set text(size:8pt)
GAE with unpriced work is a lower bound. Exhaustion is an attempted-budget
cost, not a successful DLP cost. Timeout counters/costs are missing, not zero.
Each outcome stays separate, with present-value n. Raw phase/native counters,
memory and unqualified elapsed observations remain in the native summary and
CSV. Repeated baseline and candidate targets are not pooled.

== Censored observations (one row per attempt)
#set text(size:7pt)
#table(columns:(.5fr,.5fr,1.5fr,.7fr,1fr,1fr),inset:3pt,
[Graph], [Arm], [Record id], [Final GAE], [Raw process s], [Peak RSS KiB],
[three], [bsgs], [ECR1hdd5aa78ac1bc], [missing], [3.011223532], [166736.0],
[three], [bsgs], [ECR1h1b9c3c99a3b7], [missing], [3.016101743], [166728.0],
)
#set text(size:8pt)
These process observations are not completed solve times. The final-cost table above contains present reports only; censored costs stay missing.

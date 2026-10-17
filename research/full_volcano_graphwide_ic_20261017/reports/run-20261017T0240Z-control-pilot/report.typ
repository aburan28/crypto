#set page(paper: "a4", margin: 13mm)
#set text(size: 8.3pt)
#set heading(numbering: none)
= Volcano index calculus: measured control and pilot
Run 20261017-control-pilot. Counts and correctness diagnostic; no runtime or speedup claim.

#image("relation-yield.svg", width: 100%)

*Actual scope:* 3 root vertices (conductor 1), 2 public targets each, 8 arms,
1 round, 0 warmups: 48 pilot attempts. Additional small-subgroup control:
12 attempts. This is not graph-wide completion and does not replace the
original A--E requirement. A separate resumable survey targets all 24,823
vertices and 51 strata; consult its current ledger rather than this frozen report.

*Observed:* pilot: 30 verified recoveries, 14 exhausted budgets, 4 watchdog
timeouts. Control: 8 verified, 4 exhausted. Both projected-base methods solved
both control targets. The raw-abscissa control exhausted its budget in every
pilot target. On the actual 45-bit subgroup the projected arms found zero
relations in 4,096 trials per target; no smaller subgroup was substituted.

*Validation:* all four sealed sessions pass the ecbench audit. All 38 successful
deterministic records replay identically. Supplemental native replay also
reproduces all 18 exhausted records, including every phase/counter, field
operation and base identity. Four timeouts have no final counted report;
their unknown relation counts are not treated as zero. Every recovered scalar
is checked by multiplication in the runner and auditor. The replays use the
same implementation and host; independent implementation replay is still open.

*Method:* cofactor-project each raw abscissa point, reject infinity, check
prime-subgroup membership and deduplicate negation pairs. Factor-base sizes
16 and 64; scan cap 20,000; subtraction and two-summand MITM; checked full-rank
incremental elimination. Public-target seed 20261017; order seed 1, alternate.
Cold base/table setup and verification are charged on every one-target run.
References are rho.negation and textbook BSGS, with matched rho A/A. A 3-second
child watchdog censors unfinished runs. Both methods use actual largest r.

#pagebreak()
= Counted cold costs and yields
GAE values below are *lower bounds with unpriced native work*. Each row has
two measured target attempts; min/mean/max describe those observations only.
Exhausted rows are kept separate and do not claim recovered DLP cost.

#set text(size: 6.6pt)
#table(columns: (2.6fr, 1.3fr, .9fr, .3fr, .8fr, .8fr, .8fr, .8fr, .7fr), inset: 2.5pt,
  [Source model], [Projected base/oracle], [Status], [n], [GAE min], [GAE mean], [GAE max], [Trials sum], [Relations sum],
[icv1-fp57-t536870914-6cabb0c9], [16-mitm], [verified], [2], [6,493], [6,927.0], [7,361], [5,400], [53],
[icv1-fp57-t536870914-6cabb0c9], [16-subtract], [verified], [2], [76,008], [90,295.5], [104,583], [5,400], [53],
[icv1-fp57-t536870914-6cabb0c9], [64-mitm], [verified], [2], [19,530], [19,970.0], [20,410], [2,569], [310],
[icv1-fp57-t536870914-6cabb0c9], [64-subtract], [verified], [2], [130,718], [160,668.5], [190,619], [2,569], [310],
[icv1-fp57-tm536870941-e6d0c766], [16-mitm], [exhausted], [2], [8,191], [8,204.0], [8,217], [8,192], [0],
[icv1-fp57-tm536870941-e6d0c766], [16-subtract], [exhausted], [2], [138,735], [138,748.0], [138,761], [8,192], [0],
[icv1-fp57-tm536870941-e6d0c766], [64-mitm], [exhausted], [2], [20,239], [20,252.0], [20,265], [8,192], [0],
[icv1-fp57-tm536870941-e6d0c766], [64-subtract], [exhausted], [2], [536,271], [536,284.0], [536,297], [8,192], [0],
[icv1-fp59-t1073741954-4dd0b1ef], [16-mitm], [verified], [2], [4,242], [4,419.0], [4,596], [598], [67],
[icv1-fp59-t1073741954-4dd0b1ef], [16-subtract], [verified], [2], [10,148], [12,328.0], [14,508], [598], [67],
[icv1-fp59-t1073741954-4dd0b1ef], [64-mitm], [verified], [2], [20,463], [20,467.5], [20,472], [532], [447],
[icv1-fp59-t1073741954-4dd0b1ef], [64-subtract], [verified], [2], [22,087], [22,588.0], [23,089], [532], [447],)
#set text(size: 8.3pt)

*Complete accounting:* native summary.json and statistics.csv include min,
max, arithmetic mean and explicit n by graph, conductor, arm, status and
phase; observations.jsonl retains record/workload ids and all phase-native
counters. Raw sessions retain every field operation, phase, setup, matrix
rank/dependency, target window and failure. pricing-band.json applies common
0/0.1/1 GAE native-work weights to IC and reference arms alike. No calibrated
versus uncalibrated speed ranking is inferred here.

*Host and limits:* Intel Xeon Platinum 8573C, Linux x86-64 container, four-CPU
quota and 16 GiB RAM. ecbench reserves CPU 4; every measured record earns L0:
unoptimized build, movable-thread/PSI and host-isolation limitations. Raw
in-process and end-to-end elapsed observations are retained. None is an
admitted optimized runtime distribution. Two targets/one round also do not
qualify wall-clock comparisons. PMU instructions/cycles were unavailable.

*Provenance:* projected-base source/tests and first runner commit 4eecf6b53;
exhaustion replay and survey ordering commit ee54ad603. Exact ecbench SHA-256:
#text(size: 6.7pt)[1869ca2e7be271794df8128c08a6dca0201038e2c81a3e0ce6434f7ec0d8379d].
Per-campaign manifests bind source, binary, queue and input hashes. Exact
control/pilot protocols, frozen specs, plans, records and replay receipts
accompany this report. Compiler/test exit logs are retained. New subgroup
base tests pass 2/2, including every nonzero control logarithm without a
planted-answer input. Conductor reservation was unavailable (exit 127).

*Open requirements:* exact per-vertex representation registry, optimized and
qualified full-matrix execution, all targets/rounds/base sizes, transported
edge experiments, composite-order PH, independent implementation replay and
finite-field IC. Actual embedding degrees are 30,875; 3,953,777,913,734; and
2,126. The existing supersingular degree-2 teaching MOV backend cannot
execute finite-field IC in those destinations. This is an implementation
gap, not proof that mathematical reductions do not exist. Original A/C/D/E
full generation remains storage-blocked; no original case is replaced here.

No new mathematical theorem or fastest-method claim is made. The established
IC/rho frontier remains unchanged. These results establish the measured
control/pilot behavior, not unrun strata or vertices.

#set page(paper: "a4", margin: 13mm)
#set text(size: 8.3pt)
#set heading(numbering: none)
= Volcano IC: all strata, 99-vertex checkpoint
Run 20261017T0337Z. Counted diagnostic; no admitted runtime or speedup claim.

#image("all-strata-yield.svg", width: 100%)

#table(columns: (1fr,1fr,1fr,1fr,1fr,1fr), inset: 3pt,
 [Graph key], [Audited vertices], [Planned vertices], [Audited strata], [Planned strata], [Largest End(E) conductor],
[mixed], [44], [10070], [28], [28], [1728],
[three], [25], [6561], [9], [9], [6561],
[two], [30], [8192], [14], [14], [8192],)

*Observed:* 1,584 measured attempts: 1,065 verified recoveries, 422 exhausted
budgets, 97 watchdog timeouts. All 51 conductor strata have measured vertices;
only 99 of 24,823 vertices have been audited. Stratum coverage is not complete
graph coverage. Maximal sampled conductors are 8192, 6561 and 1728.

*Integrity:* every standard audit receipt passes; successful/exhausted
deterministic reports replayed when collected. The stopped checkpoint inventory
rechecks every receipt/file hash, all 1,584 record seals, exact replay populations
and ledger/queue identities. It is a receipt-binding check, not a fresh or
independent solver replay. Exact inventory and raw frozen files accompany this report.

#pagebreak()
= Yield definitions and complete stratum tables
*Correction to earlier wording:* the earlier reports called 301,715 source
trials “attempted decompositions.” Trials include duplicate queries skipped
before the oracle. Their recorded values and historical PDFs are unchanged.
The new summary retains both denominators: 2,223,204 source trials and 2,075,145
actual decomposition calls, producing 57,176 relations across 990 IC reports.
For these guarded classic-walk arms, calls equal trials minus repeated queries
skipped. Identity skips occur before trial increment. All 594 two-summand MITM
records independently agree with their lookup counters. A zero/unknown call
denominator has no yield; timeout counters are never assumed zero.

The legacy field relation-yield remains a per-trial alias. Summary schema v2
adds explicit trial-relation-yield and oracle-call-relation-yield, numerator,
denominator, sample counts and definitions. Global totals across distinct
subgroups/bases are descriptive accounting, not a common success probability.

The tables show only size-64 projected MITM, with every conductor stratum.
GAE values are lower bounds with unpriced native work. Exhausted rows are
attempt costs, not successful DLP costs. Means use present observations and
the explicit n; multiple vertices per stratum are sampled, not a census.
Native CSV/JSON contain every arm, phase, status, range, matrix rank/dependency,
memory and raw elapsed observation. Exact generator representations remain
explicit and unregistered where noted; no new curve-name alias is inferred.

== icv1-fp57-t536870914-6cabb0c9
#set text(size: 7pt)
#table(columns: (.6fr,.9fr,.4fr,.9fr,.9fr,.9fr,1fr,1fr), inset: 2.5pt,
[Conductor], [Status], [n], [GAE min], [GAE mean], [GAE max], [Mean/trial], [Mean/call],
[1], [verified], [2], [19,530], [19,970.0], [20,410], [0.119769], [0.121079],
[2], [verified], [2], [19,505], [19,546.5], [19,588], [0.114059], [0.114919],
[4], [verified], [4], [19,466], [20,048.5], [20,714], [0.121453], [0.122819],
[8], [verified], [8], [19,467], [19,888.2], [20,857], [0.119064], [0.120233],
[16], [verified], [16], [19,210], [19,796.4], [20,524], [0.126336], [0.127545],
[32], [verified], [12], [19,434], [19,969.1], [20,954], [0.125741], [0.127051],
[64], [verified], [2], [19,797], [19,997.0], [20,197], [0.129674], [0.130887],
[128], [verified], [2], [19,324], [19,538.5], [19,753], [0.130729], [0.132345],
[256], [verified], [2], [19,260], [20,026.5], [20,793], [0.132167], [0.132942],
[512], [verified], [2], [19,514], [19,930.5], [20,347], [0.125041], [0.126247],
[1024], [verified], [2], [19,921], [20,062.0], [20,203], [0.121524], [0.122809],
[2048], [verified], [2], [19,076], [19,585.0], [20,094], [0.120734], [0.121480],
[4096], [verified], [2], [19,609], [20,067.0], [20,525], [0.120562], [0.121898],
[8192], [verified], [2], [19,369], [19,846.5], [20,324], [0.143079], [0.144412],)
#set text(size: 8.3pt)

== icv1-fp57-tm536870941-e6d0c766
#set text(size: 7pt)
#table(columns: (.6fr,.9fr,.4fr,.9fr,.9fr,.9fr,1fr,1fr), inset: 2.5pt,
[Conductor], [Status], [n], [GAE min], [GAE mean], [GAE max], [Mean/trial], [Mean/call],
[1], [exhausted], [2], [20,239], [20,252.0], [20,265], [0.000000], [0.000000],
[3], [exhausted], [4], [20,245], [20,270.0], [20,298], [0.000000], [0.000000],
[9], [exhausted], [12], [20,240], [20,273.7], [20,327], [0.000000], [0.000000],
[27], [exhausted], [22], [20,226], [20,265.6], [20,299], [0.000000], [0.000000],
[81], [exhausted], [2], [20,302], [20,319.0], [20,336], [0.000000], [0.000000],
[243], [exhausted], [2], [20,284], [20,299.5], [20,315], [0.000000], [0.000000],
[729], [exhausted], [2], [20,278], [20,283.0], [20,288], [0.000000], [0.000000],
[2187], [exhausted], [2], [20,253], [20,270.0], [20,287], [0.000000], [0.000000],
[6561], [exhausted], [2], [20,283], [20,286.0], [20,289], [0.000000], [0.000000],)
#set text(size: 8.3pt)

#pagebreak()

== icv1-fp59-t1073741954-4dd0b1ef
#set text(size: 7pt)
#table(columns: (.6fr,.9fr,.4fr,.9fr,.9fr,.9fr,1fr,1fr), inset: 2.5pt,
[Conductor], [Status], [n], [GAE min], [GAE mean], [GAE max], [Mean/trial], [Mean/call],
[1], [verified], [2], [20,463], [20,467.5], [20,472], [0.840387], [0.868095],
[2], [verified], [6], [20,012], [20,126.2], [20,241], [0.847810], [0.863518],
[3], [verified], [8], [20,024], [20,270.8], [20,595], [0.846633], [0.867332],
[4], [verified], [10], [19,899], [20,159.5], [20,527], [0.854733], [0.873263],
[6], [verified], [16], [19,900], [20,144.6], [20,550], [0.836148], [0.854720],
[8], [verified], [2], [20,050], [20,250.0], [20,450], [0.844634], [0.866823],
[9], [verified], [2], [20,080], [20,115.0], [20,150], [0.865867], [0.881503],
[12], [verified], [2], [19,813], [19,988.0], [20,163], [0.826053], [0.847380],
[16], [verified], [2], [20,142], [20,300.0], [20,458], [0.821076], [0.859763],
[18], [verified], [2], [20,008], [20,008.0], [20,008], [0.848025], [0.858225],
[24], [verified], [2], [20,086], [20,297.0], [20,508], [0.864248], [0.886644],
[27], [verified], [2], [19,951], [19,992.0], [20,033], [0.829576], [0.845611],
[32], [verified], [2], [19,880], [20,051.5], [20,223], [0.806809], [0.822622],
[36], [verified], [2], [19,875], [19,990.0], [20,105], [0.851496], [0.861641],
[48], [verified], [2], [20,241], [20,297.5], [20,354], [0.829370], [0.859703],
[54], [verified], [2], [20,100], [20,328.5], [20,557], [0.782841], [0.801320],
[64], [verified], [2], [19,996], [20,194.5], [20,393], [0.855256], [0.876868],
[72], [verified], [2], [19,964], [20,203.5], [20,443], [0.827778], [0.845433],
[96], [verified], [2], [19,919], [19,929.5], [19,940], [0.846732], [0.857399],
[108], [verified], [2], [20,273], [20,492.0], [20,711], [0.821501], [0.849152],
[144], [verified], [2], [19,922], [20,056.0], [20,190], [0.849713], [0.861567],
[192], [verified], [2], [20,115], [20,186.5], [20,258], [0.843428], [0.864101],
[216], [verified], [2], [20,022], [20,022.5], [20,023], [0.846824], [0.859916],
[288], [verified], [2], [20,118], [20,180.0], [20,242], [0.843369], [0.861071],
[432], [verified], [2], [19,933], [20,043.0], [20,153], [0.838037], [0.854658],
[576], [verified], [2], [20,153], [20,259.0], [20,365], [0.865327], [0.885298],
[864], [verified], [2], [19,900], [19,942.0], [19,984], [0.861202], [0.868679],
[1728], [verified], [2], [20,215], [20,372.5], [20,530], [0.847351], [0.868859],)
#set text(size: 8.3pt)

*Protocol and provenance:* projected bases 16/64 with MITM/subtraction,
unchanged raw-base MITM, rho.negation, textbook BSGS and rho A/A. Two public
targets per vertex, one alternating round, zero warmups, target seed 20261017,
order seed 1, 4,096 source trials and a 3-second child watchdog. Actual largest
prime subgroups are retained, including the 45-bit subgroup. Cold solver setup,
projection and verification are charged. All original measured sessions remain
byte-frozen. Source commits 4eecf6b53 and ee54ad603; complete stopped campaign,
native source, queue, ledger, receipts and frozen runner are in the archive.

*Host and limits:* Intel Xeon Platinum 8573C, Linux x86-64, four-CPU quota,
16 GiB RAM. Dev-profile ecbench hash
#text(size: 6.5pt)[1869ca2e7be271794df8128c08a6dca0201038e2c81a3e0ce6434f7ec0d8379d].
All wall observations are L0/debug and unqualified. No builds/test suites ran
during the resumed measured block; it was stopped before the release build.
Raw process time includes orchestration; no optimized speedup is inferred.

The full million-trial/eight-target/five-round/four-size matrix, original
A--E full graph request, all-vertex completion, exact representation registry,
transported-target benchmarks, composite PH, actual finite-field IC and
independent implementation replay remain open. Exact embedding degrees
30,875; 3,953,777,913,734; and 2,126 exceed the implemented MOV teaching backend.
This is an implementation gap, not proof that mathematical reductions do
not exist. A separately named release-profile survey must pass matched count
checks and audits before launch. Canonical fastest-method claims are unchanged.

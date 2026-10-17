#set page(paper:"a4",margin:13mm)
#set text(size:8.5pt)
= Repaired volcano survey checkpoint
Run 20261017T1849Z. 4,037 audited vertices, 64,592 accepted attempts.
#image("all-strata-yield.svg",width:100%)
#table(columns:(1fr,1fr,1fr,1fr,1fr),inset:3pt,
 [Graph], [Audited vertices], [Planned vertices], [Audited strata], [Planned strata],
[mixed], [1356], [10070], [28], [28],
[three], [1338], [6561], [9], [9],
[two], [1343], [8192], [14], [14],)

Observed outcomes: 41,573 verified recoveries,
20,132 exhausted attempts, 2,887 timeouts.
Native source trials 103,036,278; actual decomposition calls
98,188,909; equations 2,098,636.

These pooled counts span different subgroups/bases and are accounting totals,
not a common success probability. Multiple outcomes and all missing/censored
observations remain distinct. Stratum samples are not complete populations.

The environment returned without the previous worker. One in-flight attempt
was preserved and marked interrupted by the unchanged driver's zero-job stop;
it is excluded from accepted means. Exact termination cause was not observed.
No sealed session, binary, queue or mathematical budget was changed.

The native inventory checks record seals, receipt/file hashes, ledger/queue
identities and exact existing successful/exhausted replay coverage. This is
receipt binding, not fresh or independent solver replay. The entire stopped
campaign, including the unaccepted interruption, accompanies the inventory.

#pagebreak()
= Stratum distributions: recorded size-64 MITM
Trial yield includes repeated queries skipped before the oracle; actual-call
yield excludes them. Identity skips occur before trial increment. Zero or
unknown denominators have no yield. Timeout final counters are unknown, not
zero. GAE with unpriced work is a lower bound; exhaustion costs an attempt,
not a successful DLP. Native summary/CSV retain all other arms, phases and
present-value n. The tables below use only recorded size-64 MITM reports.

== icv1-fp57-t536870914-6cabb0c9
#set text(size:7pt)
#table(columns:(.65fr,.9fr,.45fr,.9fr,.9fr,.9fr,1fr,1fr),inset:2.5pt,
[Conductor], [Outcome], [n], [GAE min], [GAE mean], [GAE max], [Mean/trial], [Mean/call],
[1], [verified], [2], [19,530], [19,970.0], [20,410], [0.119769], [0.121079],
[2], [verified], [2], [19,505], [19,546.5], [19,588], [0.114059], [0.114919],
[4], [verified], [4], [19,466], [20,048.5], [20,714], [0.121453], [0.122819],
[8], [verified], [8], [19,467], [19,888.2], [20,857], [0.119064], [0.120233],
[16], [verified], [16], [19,210], [19,796.4], [20,524], [0.126336], [0.127545],
[32], [verified], [32], [19,063], [19,777.5], [20,954], [0.125379], [0.126577],
[64], [verified], [64], [18,937], [19,903.8], [21,779], [0.122785], [0.123939],
[128], [verified], [128], [18,907], [19,764.5], [21,045], [0.126269], [0.127420],
[256], [verified], [256], [18,833], [19,875.1], [21,476], [0.124160], [0.125334],
[512], [verified], [512], [18,763], [19,854.7], [21,859], [0.124704], [0.125875],
[1024], [verified], [1024], [18,691], [19,864.5], [22,854], [0.124112], [0.125262],
[2048], [verified], [634], [18,709], [19,847.9], [22,184], [0.124286], [0.125440],
[4096], [verified], [2], [19,609], [20,067.0], [20,525], [0.120562], [0.121898],
[8192], [verified], [2], [19,369], [19,846.5], [20,324], [0.143079], [0.144412],)
#set text(size:8.5pt)

== icv1-fp57-tm536870941-e6d0c766
#set text(size:7pt)
#table(columns:(.65fr,.9fr,.45fr,.9fr,.9fr,.9fr,1fr,1fr),inset:2.5pt,
[Conductor], [Outcome], [n], [GAE min], [GAE mean], [GAE max], [Mean/trial], [Mean/call],
[1], [exhausted], [2], [20,239], [20,252.0], [20,265], [0.000000], [0.000000],
[3], [exhausted], [4], [20,245], [20,270.0], [20,298], [0.000000], [0.000000],
[9], [exhausted], [12], [20,240], [20,273.7], [20,327], [0.000000], [0.000000],
[27], [exhausted], [36], [20,226], [20,271.4], [20,316], [0.000000], [0.000000],
[81], [exhausted], [108], [20,232], [20,278.6], [20,340], [0.000000], [0.000000],
[243], [exhausted], [324], [20,198], [20,274.7], [20,342], [0.000000], [0.000000],
[729], [exhausted], [972], [20,205], [20,275.2], [20,337], [0.000000], [0.000000],
[2187], [exhausted], [1216], [20,212], [20,274.4], [20,337], [0.000000], [0.000000],
[6561], [exhausted], [2], [20,283], [20,286.0], [20,289], [0.000000], [0.000000],)
#set text(size:8.5pt)

#pagebreak()

== icv1-fp59-t1073741954-4dd0b1ef
#set text(size:7pt)
#table(columns:(.65fr,.9fr,.45fr,.9fr,.9fr,.9fr,1fr,1fr),inset:2.5pt,
[Conductor], [Outcome], [n], [GAE min], [GAE mean], [GAE max], [Mean/trial], [Mean/call],
[1], [verified], [2], [20,463], [20,467.5], [20,472], [0.840387], [0.868095],
[2], [verified], [6], [20,012], [20,126.2], [20,241], [0.847810], [0.863518],
[3], [verified], [8], [20,024], [20,270.8], [20,595], [0.846633], [0.867332],
[4], [verified], [12], [19,899], [20,210.8], [20,534], [0.852565], [0.872370],
[6], [verified], [24], [19,899], [20,142.7], [20,550], [0.835873], [0.854372],
[8], [verified], [24], [19,807], [20,142.4], [20,648], [0.835348], [0.854639],
[9], [verified], [24], [19,842], [20,207.2], [20,878], [0.847141], [0.867336],
[12], [verified], [48], [19,783], [20,183.1], [20,749], [0.842426], [0.861242],
[16], [verified], [48], [19,796], [20,210.0], [20,727], [0.834727], [0.855684],
[18], [verified], [72], [19,863], [20,202.9], [20,836], [0.837101], [0.859458],
[24], [verified], [96], [19,835], [20,193.0], [21,445], [0.846245], [0.865569],
[27], [verified], [72], [19,786], [20,146.4], [20,934], [0.844967], [0.864468],
[32], [verified], [96], [19,832], [20,191.3], [20,771], [0.846372], [0.865296],
[36], [verified], [144], [19,723], [20,156.5], [20,928], [0.842766], [0.860938],
[48], [verified], [192], [19,734], [20,208.9], [21,208], [0.843601], [0.863686],
[54], [verified], [216], [19,757], [20,230.8], [21,129], [0.838622], [0.858583],
[64], [verified], [92], [19,829], [20,149.4], [21,106], [0.845215], [0.864308],
[72], [verified], [288], [19,741], [20,178.8], [21,201], [0.847042], [0.866266],
[96], [verified], [180], [19,776], [20,215.6], [20,971], [0.841509], [0.861867],
[108], [verified], [432], [19,785], [20,225.3], [21,246], [0.840016], [0.861047],
[144], [verified], [264], [19,829], [20,211.4], [21,051], [0.844375], [0.865404],
[192], [verified], [2], [20,115], [20,186.5], [20,258], [0.843428], [0.864101],
[216], [verified], [360], [19,797], [20,206.2], [21,109], [0.844457], [0.865507],
[288], [verified], [2], [20,118], [20,180.0], [20,242], [0.843369], [0.861071],
[432], [verified], [2], [19,933], [20,043.0], [20,153], [0.838037], [0.854658],
[576], [verified], [2], [20,153], [20,259.0], [20,365], [0.865327], [0.885298],
[864], [verified], [2], [19,900], [19,942.0], [19,984], [0.861202], [0.868679],
[1728], [verified], [2], [20,215], [20,372.5], [20,530], [0.847351], [0.868859],)
#set text(size:8.5pt)

Protocol: actual largest prime subgroups; two public targets/model, one
alternating round, zero warmups, bases 16/64, raw-base control, MITM/subtraction,
rho.negation, BSGS and rho A/A. 4,096 source trials, 3-second child watchdog.
Cold setup/projection/verification charged; unpriced work remains explicit.
The full eight-target/five-round/million-trial matrix remains unexecuted.

Frozen release runner 35d89073f4d666b094c213a1d31b0c819b63600123e067c111af3c947d7609a9. Source 2bf7700b7, unchanged
queue and driver. Earlier dev99, stopped release70, control/pilot and original
C finite-field diagnostic are separate datasets; repeated workloads and
different binaries are never pooled here. No performance improvement follows
from this checkpoint's larger sample count.

Wall observations are unqualified and descriptive. Every capsule retains its
host/environment observation; any resumed runtime epoch must remain identifiable.
No admitted runtime, external fastest-method, scaling or sub-rho claim is made.
Exact representation registration, transported targets, composite PH,
independent arithmetic replay, original A--E full graphs and large-extension
finite-field backends remain open. Original C's r=7 k=1 diagnostic is separate.
Complete graph topology is for the declared fields and finite prime degrees,
not every other prime or infinitely many compositions.


== This checkpoint's clean pause
The current invocation stopped normally at a vertex boundary on an explicit
STOP marker, at 18:48:59 UTC on 17 October. It completed 1,102 additional
vertices after the earlier 2,935 checkpoint, yielding 4,037 accepted vertices.
The earlier interruption is historical; its original attempt remains retained
and attempt 2 completed. This is not a new failed or interrupted capsule.

The exact queue, driver and solver binary are unchanged. Method search budgets
remain 4,096 trials and a three-second child watchdog. This inventory spans
separate runtime epochs: aggregate operation/yield summaries are descriptive;
all clock fields and their aggregate means remain unqualified and mixed-epoch.
Do not use them for admitted timing comparisons or pool them with the older
dev99/release70 or the separate original-C full-point-order diagnostic.

A later resume has its own receipt and outer resource guards. No future run,
S3 upload or full-vertex completion is implied by this frozen snapshot.

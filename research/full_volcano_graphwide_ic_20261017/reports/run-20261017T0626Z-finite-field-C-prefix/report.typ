#set page(paper: "a4", margin: 13mm)
#set text(size: 9pt)
= Original C: finite-field reduction diagnostic
Source model icv1-fp61-t2-fefb1f9c. Public synthetic inputs, 17 October 2026.
This report is descriptive; no admitted runtime or scaling claim.

The original C model has p=1637958981148213249,
a=571607911220311176, b=29329524900995129 and
N=1637958981148213248. Native complete factorization gives
N=2^26 × 3^20 × 7. Its largest prime subgroup is 7 and embedding degree is 1.
The 61-bit field is not a 61-bit prime-order DLP.

#image("route.svg", width: 100%)

The shifted-divisor Miller implementation computes standard reduced Tate
evaluations. Finite tests check all six nonzero control/C scalar inputs and
all nontrivial destination generators/targets; these are implementation tests,
not a new bilinearity theorem or independent implementation proof. Recovered
factor-base logs, the destination answer and the original EC equation are
checked. The solver receives neither planted scalars nor factor-base logs.

== Original-input arithmetic inventory (not five pairing runs)
#table(columns: (.4fr,1.4fr,1.6fr,1.8fr), inset: 4pt,
 [Case], [Largest prime r], [Exact embedding degree], [Dense u64 bytes / field element],
[A], [716481758899867], [238827252966622], [1910618023732976],
[B], [15024241], [341460], [2731680],
[C], [7], [1], [8],
[D], [101139959083], [101139959082], [809119672656],
[E], [87131751841], [17426350368], [139410802944],
)
Dense coefficient sizes describe that specific layout only. They are not
universal lower bounds, impossibility proofs or estimates for other algorithms.
The A/B/D/E finite-field routes have not been executed by this diagnostic.

== Coverage and retained limitations
#table(columns: (2fr,1fr,1fr,1fr), inset: 4pt,
 [Population], [Requested targets], [Verified], [Failed],
[Original source model], [2], [2], [0],
[Stored prefix], [3738], [3738], [0],
)
The prefix contains 1,869 stored vertices, 1,024 expanded vertices and 4,903
replayed maps. Its frontier is not exhausted; global stratum populations remain
NOT-ESTABLISHED. The separate conservative local certificate is described on
page 4. Original C's expected class has 1,934,901,248 vertices.
The source-model and prefix populations stay separate, including any repeated
target. Two public searches may select the same subgroup point; they are not
independent random samples. Other prime-power components are not measured.

#pagebreak()
= Recorded prefix work and elapsed observations
#image("prefix-work.svg", width: 100%)

Relation yield has two explicit meanings: number of collected equations per
destination query (several per query are allowed), and fraction of queries
with any relation. Identity queries are skipped separately. Relations satisfy
the checked destination equation before matrix insertion. Full rank is three.
#table(columns: (2fr,.5fr,1fr,1fr,1fr), inset: 4pt,
 [Metric], [n], [Min], [Mean], [Max],
[query-hit-yield], [3738], [1.0], [1.000000], [1.0],
[relations-per-query], [3738], [1.0], [1.708266], [2.5],
)

Native per-target end-to-end elapsed interval: setup/order factorization,
point search/projection, pairing, cold factor base, collection/matrix and
destination/EC verification. Input file loading, batch serialization and
output I/O are outside that per-target interval. Bash process elapsed/user/system
logs are separate; peak RSS was not measured because GNU time is absent. All observations are unqualified; there is no external
reference or fastest-method claim.
#table(columns: (.5fr,1fr,1fr,1fr), inset: 4pt,
 [n], [Min ns], [Mean ns], [Max ns],
[3738], [121110], [211175.290], [3620010],
)

The archive retains both executions, exact binary, source, protocol, fixed
prefix DB and input hashes, native summaries and per-target data. The native
comparison removes only elapsed time fields, then compares every remaining
value, including source/model identities, outcomes, relations and operation
ledgers. It is same-binary reproduction, not independent implementation replay.
See root-replay-match.json and prefix-replay-match.json for actual checks.

Unpriced work stays explicit: order factorization internals; square-root
internals (calls counted); allocation/comparisons; matrix inversions (matrix
entry operations counted). Explicit Fp inversions are expanded into their
multiplication/square counts: pricing those inversions again would double count.
Native EC field counts are separate from the explicit Miller/destination work.
There is no common-cost conversion or complete end-to-end GAE claim.

The three complete additional classes are separate datasets. Full A--E graph
completion, all-vertex/full-matrix benchmarks, all requested strata, exact
representation registration and large-extension finite-field backends remain open. The original graph explorer
and fastest-method scoreboards are unchanged: this diagnostic establishes no
new complete graph, general algebraic identity or admitted speed comparison.

#pagebreak()
= Per-phase operation distributions (native summary)
#set text(size: 6.7pt)
#table(columns: (1fr,2.25fr,.5fr,.8fr,.9fr,.8fr), inset: 2pt,
 [Phase], [Primitive ledger], [n], [Min], [Mean], [Max],
[EC check], [ec-field/invs], [3738], [0], [1.849], [3],
[EC check], [ec-field/muls], [3738], [0], [7.716], [12],
[EC check], [ec-field/sqrs], [3738], [0], [3.188], [5],
[EC check], [ec-group/adds], [3738], [0], [0.509], [1],
[EC check], [ec-group/doubles], [3738], [0], [1.339], [2],
[EC check], [ec-group/scalar-mults], [3738], [1], [1.000], [1],
[EC check], [explicit-fp/additions-subtractions], [3738], [0], [0.000], [0],
[EC check], [explicit-fp/exponentiations], [3738], [0], [0.000], [0],
[EC check], [explicit-fp/inversions], [3738], [0], [0.000], [0],
[EC check], [explicit-fp/multiplications], [3738], [0], [0.000], [0],
[EC check], [explicit-fp/squarings], [3738], [0], [0.000], [0],
[EC check], [miller-steps], [3738], [0], [0.000], [0],
[EC check], [shift-candidates], [3738], [0], [0.000], [0],
[EC check], [sqrt-calls-unexpanded], [3738], [0], [0.000], [0],
[FF solve/checks], [ec-field/invs], [3738], [0], [0.000], [0],
[FF solve/checks], [ec-field/muls], [3738], [0], [0.000], [0],
[FF solve/checks], [ec-field/sqrs], [3738], [0], [0.000], [0],
[FF solve/checks], [ec-group/adds], [3738], [0], [0.000], [0],
[FF solve/checks], [ec-group/doubles], [3738], [0], [0.000], [0],
[FF solve/checks], [ec-group/scalar-mults], [3738], [0], [0.000], [0],
[FF solve/checks], [explicit-fp/additions-subtractions], [3738], [0], [0.000], [0],
[FF solve/checks], [explicit-fp/exponentiations], [3738], [14], [14.000], [14],
[FF solve/checks], [explicit-fp/inversions], [3738], [5], [5.000], [5],
[FF solve/checks], [explicit-fp/multiplications], [3738], [322], [334.924], [347],
[FF solve/checks], [explicit-fp/squarings], [3738], [533], [535.770], [538],
[FF solve/checks], [miller-steps], [3738], [0], [0.000], [0],
[FF solve/checks], [shift-candidates], [3738], [0], [0.000], [0],
[FF solve/checks], [sqrt-calls-unexpanded], [3738], [0], [0.000], [0],
[Pairing], [ec-field/invs], [3738], [14], [14.000], [14],
[Pairing], [ec-field/muls], [3738], [52], [52.000], [52],
[Pairing], [ec-field/sqrs], [3738], [22], [22.000], [22],
[Pairing], [ec-group/adds], [3738], [6], [6.000], [6],
[Pairing], [ec-group/doubles], [3738], [8], [8.000], [8],
[Pairing], [ec-group/scalar-mults], [3738], [0], [0.000], [0],
[Pairing], [explicit-fp/additions-subtractions], [3738], [80], [84.195], [120],
[Pairing], [explicit-fp/exponentiations], [3738], [28], [28.000], [28],
[Pairing], [explicit-fp/inversions], [3738], [26], [26.000], [26],
[Pairing], [explicit-fp/multiplications], [3738], [1270], [1274.195], [1310],
[Pairing], [explicit-fp/squarings], [3738], [1692], [1694.097], [1712],
[Pairing], [miller-steps], [3738], [8], [8.000], [8],
[Pairing], [shift-candidates], [3738], [2], [4.097], [22],
[Pairing], [sqrt-calls-unexpanded], [3738], [2], [4.097], [22],
[Setup], [ec-field/invs], [3738], [152], [174.292], [511],
[Setup], [ec-field/muls], [3738], [658], [754.127], [2213],
[Setup], [ec-field/sqrs], [3738], [270], [309.473], [908],
[Setup], [ec-group/adds], [3738], [36], [41.111], [116],
[Setup], [ec-group/doubles], [3738], [118], [136.207], [403],
[Setup], [ec-group/scalar-mults], [3738], [4], [4.319], [9],
[Setup], [explicit-fp/additions-subtractions], [3738], [4], [9.331], [42],
[Setup], [explicit-fp/exponentiations], [3738], [0], [0.000], [0],
[Setup], [explicit-fp/inversions], [3738], [0], [0.000], [0],
[Setup], [explicit-fp/multiplications], [3738], [4], [9.331], [42],
[Setup], [explicit-fp/squarings], [3738], [2], [4.665], [21],
[Setup], [miller-steps], [3738], [0], [0.000], [0],
[Setup], [shift-candidates], [3738], [0], [0.000], [0],
[Setup], [sqrt-calls-unexpanded], [3738], [2], [4.665], [21],
)
#set text(size: 8pt)
Source: #raw("execution/finite-field-k1-evidence"),
#raw("execution/finite-field-k1-summary.json"),
#raw("execution/FINITE_FIELD_K1_PROTOCOL.md") and the delivery manifest. Counts cover
present verified observations only; failed rows remain in raw data. All raw
target and primitive fields are also exported in prefix-observations.csv.

#pagebreak()
= Local conductor bounds: partial coverage of original C
#image("local-strata.svg", width: 100%)

The local certificate labels 517 vertices across
24 conductor strata, from
2,985,984 to
967,458,816. The source model's exact local
conductor is 967,458,816. The other
1352 vertices retain unknown conductor;
2704 target records remain in that separate bucket.

This uses the existing ordinary-volcano floor-distance correspondence. Complete
kernel multiplicities identify known floors. Distance to known floors bounds
height above; distance to known floors or unexpanded frontier bounds it below.
Only equality of finite bounds produces an exact local height/conductor.
Unexpanded nodes are possible floors for the bound, never asserted floors.
The frozen prefix's global stratum receipt and database are unchanged.

The second labeling implementation passed all checks against the three
complete graphs at visibility cutoffs 0, 1, 17, 1024 and all. These are finite
calibrations using the same standard volcano theorem and isogeny arithmetic,
not a formal proof or independent arithmetic implementation. The full-class
population and all 165 strata remain unestablished. Local annotation is study
preparation and is not used by the DLP solver or charged to its elapsed interval.

Native per-stratum summaries keep all phases, statuses, min/max/arithmetic
means and present-value n. An unknown conductor is a null bucket, not conductor
zero or an omitted observation. No causal effect of conductor on cost is
inferred from these differing points and sample sizes.

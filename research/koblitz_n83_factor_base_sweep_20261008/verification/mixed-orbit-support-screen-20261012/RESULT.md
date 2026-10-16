# Exact support screen after the first mixed-orbit search

This is **integer counting from already audited point counts**, with no new
retained-N83 solver run and no additional experimental charge. The [import
receipt](../mixed-orbit-search-15000-20261012/evidence-01/run/verify.stdout)
fixes 196,212 small-base points and 143,756 disjoint residual points in the
K=1,182/K=2,048 nested pair. The [search audit](../mixed-orbit-search-15000-20261012/audit.json)
binds those counts to source `f0a1e55f7cb6d965b0791ef68cb636ec8b8a3a1a`
and to the one-trial `UNKNOWN_search_cap` result. The [screen](screen.json)
records the exact receipt hashes and six counts; [support_screen.py](support_screen.py)
recomputes them without opening an S3 object or invoking a solver.

For `m` total summands and exactly `j` residual summands, allowing repeated
points within each disjoint domain, the number of unordered point multisets is
`C(B+m-j-1,m-j) × C(L+j-1,j)`. Every multiset has one group sum. On a uniformly
chosen nonidentity target in the primary subgroup of order
`2417851639230796216685689`, its expected decomposition multiplicity is **at
most** the count divided by `r-1`: multisets summing to the identity contribute
zero. Markov clips the corresponding hit-probability ceiling at one. Neither
bound estimates solver success on the frozen public target.

| Split, small + residual | Exact multisets | Count / (`r-1`) | Hit-probability ceiling |
| --- | ---: | ---: | ---: |
| 5 + 0 comparator | 2,423,643,343,511,761,316,264,808 | 1.00239539 | 1 |
| 4 + 1 | 8,878,309,426,598,156,108,089,140 | 3.67198271 | 1 |
| 3 + 2, measured SAT shape | 13,009,394,065,076,279,872,900,344 | 5.38055928 | 1 |
| 6 + 0 comparator | 79,260,004,322,307,878,365,588,638,556 | 32,781.1695 | 1 |
| 5 + 1 | 348,413,272,489,876,759,780,963,738,848 | 144,100.352 | 1 |
| 4 + 2 | 638,159,564,119,735,563,815,285,249,490 | 263,936.61 | 1 |

The 5+0 value independently agrees with the earlier K=1,182 nominal m=5
capacity in [SOLVER_GATES.md](../../SOLVER_GATES.md). All six Markov ceilings
are vacuous. Thus the point counts neither exclude a decomposition nor explain
the zero models at 5,000 conflicts. The empirical run resolved the model
construction and memory gates but censored SAT search before a model. The
source's first identity mask is the generic all-affine case; later masks
describe intermediate sums equal to the identity. Dividing a fixed conflict
budget evenly among those exceptional masks has no support-based rationale.

The next useful gate is solver access to an exact, independently verified N83
witness under the same finite domains, explicitly labelled as a known-answer
control. Only after that control can a new natural-target search test yield and
graph closure. Any timing comparison of bases still needs full cold setup,
natural rank, sparse linear algebra, a verified individual logarithm, and a
matched rho run. This screen supplies none of those phases and selects no base.

Validation: fresh `support_screen.py` output matched `screen.json` byte for
byte; the prior retained-run audit passed again; the study and boundary Python
suites passed 46 and 43 tests. The required release library suite passed with
4,000 tests and 113 ignored using four test threads. An earlier, more
contended default-concurrency run had 3,998 passes and two unrelated timed F4
tests skip too many cases; each passed alone before the full four-thread pass.

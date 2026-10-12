# Decision: K8 wins cold whole-solve instructions; online base choice stays open

The untouched 16-point confirmation selected K8 under the
[preregistered](PROTOCOL.md) **cold whole-solve Callgrind instruction** rule.
All 64 first-measured-round profiles recovered the same public-point
logarithms as the sealed native session. K16/K8 is **1.19026** with a
target-block bootstrap 95% interval **[1.15842, 1.22175]**, wholly above the
1.05 selection threshold. K8 therefore uses about 16.0% fewer instructions
than K16 for this complete `methods::solve` implementation on these points.
The duplicate K16 control's largest paired deviation is **0.00415%**, below
the frozen 2% A/A limit.

| One-target arm | Exact candidate or method ID | Actual usable base points | Folded columns | Verified profiles | Mean solve Ir | Mean Ir / √r |
|:--|:--|--:|--:|--:|--:|--:|
| Strong signed-Frobenius rho | `ECM1h1d9961ee4601` | — | — | 16/16 | 9,627,417 | 633.98 |
| IC K8 | `IC1N37Ckb0fb592PDP3mitmfrobeniuscountedRCsampleLAgaussTDpdpISO0h0eb4fc2e0d54` | 592 | 8 | 16/16 | 35,088,008 | 2,310.61 |
| IC K16 | `IC1N37Ckb0fb1184PDP3mitmfrobeniuscountedRCsampleLAgaussTDpdpISO0h44f5af6dc772` | 1,184 | 16 | 16/16 | 41,763,838 | 2,750.22 |
| Identical K16 control | same K16 candidate | 1,184 | 16 | 16/16 | 41,763,767 | 2,750.22 |

The unit throughout that table is **Callgrind-simulated user-space
instructions inside `methods::solve`**, from its pre-solve reset through its
post-solve reset. `Ir / √r` is labelled `crypto.S.callgrind_ir`, with
`r = 230603167`; it is not a group-addition equivalent and cannot be
compared with the generic-group floor. The 16 target-paired ratios of sums
are K16/K8 = **1.19026** [1.15842, 1.22175], K8/rho = **3.64459**
[3.23534, 4.13602], K16/rho = **4.33801** [3.84523, 4.95096], and
K16/K16-control = **1.00000172** [0.99999279, 1.00001077]. Intervals are
fixed-seed, 20,000-sample target-block bootstrap percentiles. Every raw
per-target count, including the A/A pairs, is in [DECISION.json](DECISION.json)
and [CENSUS.json](CENSUS.json). Rho remains cheaper in this unit; no IC/rho
crossover is established.

The complete [sealed ecbench session](sessions/mac_arm64_l0_01) has 384/384
verified executions, including 320 measured executions. Its local
[full audit](AUDIT.json) and independent Linux
[receipt](independent_validation/RECEIPT.json) reproduced all 320 measured
records exactly. The latter was produced by
[CI run 37191831916](https://github.com/aburan28/crypto/actions/runs/37191831916)
on a different environment class. The [32 one-target IC1 claim
records](candidate_claims) map both candidates to all 16 exact canonical
workloads and distinct run IDs; each pairs with strong rho on the same
public point, verifies both scalars, and passes its `vs_rho` correctness
check. Every claim remains **descriptive only** because both measured arms
earned L0, below the required L2 isolation level. The **incomplete** cold counted-GAE
comparison over five rounds was K16/K8 = 1.0175 [0.973, 1.059], so it could
not distinguish the two bases. This profile resolves that specific native-work
gap in a separate instruction unit. The two units must not be combined into
a speedup ratio.

The mechanism is consistent with a setup-versus-target tradeoff, not yet a
measured phase-level Ir attribution. K8 constructs 592 usable points and
2,344 pair-table entries per run; K16 constructs 1,184 points and 9,412
entries. In the sealed session, K8 averaged 5.6625 target attempts and
3,384.60 **counted** target GAE, versus K16's 1.1125 attempts and 477.33
counted target GAE. Thus the cold instruction result selects K8 for the
**cold preparation route**, while K16 remains a live candidate for the
primary one-target **online** question after reusable preparation. An
isolated target-only instruction/wall comparison is needed before discarding
either base at n41 or n53. The Mac session earned L0 on every run, so its wall
times are exploratory. This result supplies no isolated online wall speedup,
fully charged cold wall speedup, or n41/n53/n131 transfer claim.

The profiles used one release binary, SHA-256
`7e90fcf5520355ffd1d1352e0ce172127a803908947c73db8471372c0313a5a5`,
under Valgrind 3.22.0 on an Ubuntu 24.04 x86-64 hosted VM reporting AMD
EPYC 7763. The exact dependency lockfile is inside the raw archive, SHA-256
`53ad716677068cfc71a6e70d726a3be930c27e5314d1503f8fd9af8f9d1bb0fb`.
The [raw profile ZIP](RAW-CALLGRIND.zip) is 5,578,846 bytes, SHA-256
`4e0d602d7770deb9c7216137778c11ef0614cb2330c9dd95f74003a47ceb3000`.
Its [compact census](CENSUS.json), [host/build provenance](PROVENANCE.txt),
independent [receipt](independent_validation/RECEIPT.json), and
[decision](DECISION.json) have SHA-256 values
`3e04bd86ee012550870676737a9f65e1f1856e2c7b69d45f9ee43dd690ccf182`,
`8983f4c6c1a7689e9a44cc70d696a1e619a640c0a98a75cce355938fc4f7757e`,
`8a5f8683c8b3aea4c1e8343242d4d6319f4938fb8f6ae4a342c9c4ad3dad5baa`,
and `c4a33ffd9045a479ffe2f2a0ae2ecb60519cb293ac326d4ac84b4014f7980127`
respectively. The [native analyzer](../../examples/ecbench_k8_k16_analyze.rs)
re-hashes every raw file, re-parses every Callgrind part, checks the exact
method/workload/scalar against the sealed session, verifies the other-host
receipt, and recomputes the decision. Extracting the committed ZIP and
running that analyzer reproduces `DECISION.json` byte for byte. The profiler's
CPU feature dispatch under Valgrind was not independently pinned, so these
are simulated execution counts, not native retired instructions.

The first CI attempt
[failed before profiling](https://github.com/aburan28/crypto/actions/runs/37191621092)
because `Cargo.lock` is not tracked and the initial workflow passed
`--locked`. The corrected workflow archived its generated lockfile. No
measurement row was omitted or replaced by that build failure.

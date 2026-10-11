# Large-field population validation

The pre-run protocol fixes ten fields, 512 independent source fixtures, four GP
workers, and separate admission, direct-model, torsion, and class-witness labels.

| Check | Verified result and artifact |
| --- | --- |
| Native process completion | 512/512 population cells, 10/10 pilot cells, 10/10 large positive controls, and 68/68 small-field cells complete; every raw stderr is empty |
| Population result | 323/512 ordinary sources admitted; all 323 have a full-4 representative; 189 theorem exclusions; 0 direct weak source models |
| Class searches | 32,484 tested vertices; 54,216 degree-2 and 38,093 degree-3 evaluated search edges; 219 component closures and 104 vertex caps; admitted class labels unresolved |
| Positive witnesses | Ten large controls and four small controls pass fresh endpoint counts, exact norms, source/target point checks, and saved-route reconstruction |
| Small exact oracle | 64 independent sources belong to 42 weak classes and 22 zero classes; 38 search witnesses; seven zero and four positive classes among 11 restricted closures; all 15 exclusions are exact zeros |
| Raw/hash and identity replay | Native checker rebuilds the source freeze, raw hashes, records, and summaries from every completed receipt; 868 distinct population model identities |
| Independent model reconstruction | `model_replay.txt`: 1,326 model records and 72 saved route edges, reconstructed from exact moduli and roots; all short coefficients, j-invariants, rational halves, norms, image points, and order replays pass |
| Algebra receipts | Four identity records pass 32 recorded plus 64 extra evaluations; changed duplication polynomial is rejected by the native test |
| Norm-locus enumeration | 192,317 nonexcluded parameters across three fields; each event has exactly N−1 parameters and every union satisfies the density bound |
| One-edge theorem enumeration | All 5,264 ordered pairs in eight fields pass; 1,336 admitted cases comprise 592 already-full-4 and 744 requiring an edge; all quotient orders independently agree |
| Exact CM support preflight | Fundamental discriminant and f_pi=8 obtained; PARI `polclass` reports `overflow in t_INT-->long assignment`; receipt labels the class unresolved despite GP exit code zero |
| Native tests | Focused `report_check_runtime` bins pass, including retained populations, earlier theorems, orbit records, identity mutation controls, and documentary checks |
| Publication | Complete report source, SVG plot and diagrams, compiled PDF, canonical dashboard source and embedded copy, curve source manifest and registry |
| Broader published family | 100 constructed quadratic-h models; 46 ordinary traces beyond the cubic predicate; three complete all-x counts, explicit coordinate maps, and a frozen-input replay all pass |
| Required root library test | Exit 101: this isolated sparse checkout lacks `src/lib.rs`; raw stderr retained. Focused native checks pass independently. The PR's broader Rust and dashboard checks also have pre-existing failures and remain merge gates. |

The 512 primary sources are uniform ordered distinct nonzero root pairs over
the recorded fields. Class-uniform small-prime census precision, directly
constructed weak-curve controls, and this curve-weighted admission measurement
are distinct populations. The admitted large-field labels have no exact zero
or positive class decisions, giving identified precision interval [0,1].

These labels refer to the cubic norm-one branch. The prior-work comparison
adds a quadratic representative at trace 610 over F_(7^6), plus two depth-one
quartic controls. This corrects the interpretation of a cubic zero as a zero
for every Joux–Vitse form. The executed protocol and raw cubic data remain
frozen; current documentary labels explicitly name the cubic family.

All arithmetic and experiment orchestration use native Rust and compiled PARI.
Ruby updates documentary dashboard context only. The native presenter registers
the model identities without changing older model values. Registering these
curves exposed a missing object delimiter in the existing registry; the two
missing delimiter lines are repaired and the resulting JSON is parsed again.

Resource receipts use an Apple M4 Pro, 14 logical CPUs, 48 GiB memory,
Darwin 25.6.0 arm64, PARI/GP 2.17.3, and Rust 1.93.1. Shared-host wall
times retain their diagnostic status. The primary process wall receipt is
494,663 ms; it includes process starts, curve generation, point counting,
torsion checks, searches and verification for all 512 sources.

The run driver records GP and native executable hashes and refuses to overwrite
an earlier run. `source_freeze.txt` binds the protocol and executed sources.
The separate final evidence manifest binds the replay, presentation and report
artifacts. Historical source/evidence seals remain unchanged.

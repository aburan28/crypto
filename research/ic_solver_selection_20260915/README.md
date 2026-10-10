# Fixed IC solver selection: implemented, no routing gain in the initial portfolio

## Restore the research evidence

Large text artifacts are stored as readable chunks so each publication request
fits the automatic review limit. Restore their exact bytes before training on
the saved panels or rerunning an audit:

```sh
python3 research/ic_solver_selection_20260915/restore_evidence.py
```

The script verifies every artifact's SHA-256 hash and refuses to replace a
changed file. Use `--verify` to validate chunks without writing files.

The opt-in selector, matched query collector, cost-sensitive trainer and JSON
model integration are implemented. **The fitted model chooses pairs.** Every
measured natural query favored the prepared pair-table policy over SAT with
pair fallback. The hindsight oracle and the best constant policy have identical
cost and rank gain, so this portfolio supplies zero measured routing headroom.
The existing default remains pairs. No performance winner is promoted.

The initial interface covers per-query solver and budget selection for an
already chosen factor base. It does not select a new base or add a Gröbner
solver. The [usage guide](../../docs/ic/SOLVER_SELECTION.md) gives runnable
collection, training and inference commands.

## Contract and boundary

[contract.json](contract.json) was saved before implementation, alongside a
snapshot of the current baseline including local changes. Candidates use the
same curve, subgroup, factor base, target, summand count and prior matrix within
each panel. The policies are pairs with the full budget, SAT with 10% plus a
90% pair fallback, and SAT with 50% plus a 50% pair fallback. The native budget
is 0.25 seconds, with the existing separate SAT child watchdog. Policies retain
every failure, timeout, verification and rank test.

With B signed base points and m summands, uniform nonzero-target decomposition
probability is at most `min(1, C(B+m-1,m)/(r-1))`. Routing changes neither B nor m.
For the existing degree-131, three-summand profile, the inherited counting
ceiling remains approximately `1.209546e-29`. ML does not supply missing natural
relations. The reference is a matched rho workload in a common calibrated unit;
that unit has not been calibrated for this Python adapter, so rho/floor ratios
and S remain null.

The parent [tournament contract](parent-OPERATIONS.md)
admits Rust CPU implementations and counts full amd64 Valgrind instructions.
This Python portfolio has a new functional/evidence adapter and is **not
admitted to that performance tournament**. Its independently verified native
wall-time results are diagnostics. The existing Rust tournament's batch16
incumbent is unchanged. Future performance admission still requires the parent
20% full-instruction gate, paired 95% upper ratio below one and no cell more
than 10% worse.

## Matched query evidence

Frozen [panel plan](panel-plan.json): 288 natural queries and 16 separately
tagged planted controls, each evaluated by all three policies, for 912 policy
outcomes. Queries are uniform nonzero scalar multiples of a generator. Scalars
are not selector inputs. Each record retains the exact target, prior verified
relations, features, witness, transported row, rank gain and elapsed cost.

| Natural panel | Queries | Verified relations per policy | Rank gains per policy |
|---|---:|---:|---:|
| Degree 9, two summands | 48 | 44 | 44 |
| Degree 9, three summands | 48 | 48 | 48 |
| Degree 23, two summands | 96 | 3 | 3 |
| Degree 23, three summands | 96 | 90 | 89 |
| Total | 288 | 185 | 184 |

This is the useful success signal missing from the earlier full-size,
near-zero-yield dataset. The two-summand panel provides failures; the
three-summand panel provides frequent decompositions and an observed dependent
relation. Neither planted success nor an unresolved timeout is treated as
natural decomposition probability.

Whole signed Frobenius orbits determine the split, independently of generator
and base. Matrix histories also stay within their split. There are 195 training
queries from 116 distinct orbits, 50 validation queries and 43 confirmation
queries. The tiny degree-9 group has no confirmation orbit under the frozen
hash rule; all 43 confirmation queries are degree 23. Repeated members of a
training orbit do not satisfy the tree's minimum distinct-orbit leaf count.

The depth-three tree optimizes elapsed cost plus a non-gain penalty equal to
the best constant training cost per rank gain. It needs a 5% validation-loss
improvement to replace the constant. Observed improvement was **0%**. The
locked [model](model.json) therefore contains a single `pairs` leaf. The
[model report](model.report.json) preserves static policies, selector and
hindsight-oracle results for every split, and the fitting time. The selector and hindsight columns replay recorded policy
costs; they exclude model-inference overhead. Actual inference overhead is
charged by the separate complete-workload adapter.

Costs per query include feature extraction, policy calls, witness checking and
rank labeling. Pair construction and symbolic circuit setup are shared once
per panel and explicitly reported; these warm costs do not establish cold
speedup. Acquisition, failed collection, rejected labeling work and training
costs remain in the artifacts and are not amortized away in a claimed win.

## Independent audit and retained rejection

[audit.py](audit.py) computes the normal-to-polynomial basis map using an
independent direct coordinate product and uses the parent tournament's frozen
[polynomial-basis checker](oracle_pb.py). It checks group membership, signed
Frobenius transport, each decomposition, cancellation exclusions, prior matrix
rank, useful-rank labels, feature values and orbit assignment. A separate
pair/triple-sum oracle checks every final no-decomposition certificate.

The first audit rejected a rank label. The benchmark loop was reading stale
row values after reduction replaced its list. This did not alter the solver
witnesses, but it invalidated the labels. The original `candidate/`, `panels/`
and [rejection](audit-rejection-01.json) are retained. The corrected
`candidate-v2/` reran the **same frozen panel plan**, and a regression test now
covers dependent rows after elimination. No model was trained on the rejected
labels. The accepted [query audit](query-audit-v2.json) verifies all 304 records
and 603 successful policy witnesses, including controls.

## Complete-workload adapter and correctness

The baseline [A/A check](aa.json) completed before comparison and independently
verified identical relation transcripts. The [full plan](full-plan.json) fixes
12 cells: degrees 5, 9 and 23, two/three summands, development and fresh
confirmation seeds, two targets each, and three repetitions per variant.
Every worker starts with an empty durable campaign, a 100,000 pair-construction
budget, 512 attempts per stream/target and the same 0.25-second query budget.

[worker.py](worker.py) imports frozen source. [compare.py](compare.py) locks the
model hash, measures the whole process, preserves raw output and timeouts, and
independently checks every relation, base logarithm and returned target scalar.
It charges startup, cold precomputation, model loading/inference, failed
queries, verification, matrix work, commits, certificate export, reporting and
exit. Independent external audit time is outside that worker scope. All **72/72 processes completed**, with **144/144 target results** and
**420 relations** independently verified. The [final audit](final-audit.json)
also checked every full-workload rank label and reproduced the locked model
byte for byte through the training CLI. The
[summary](full-comparison/summary.json) and [raw comparisons](full-comparison/raw.jsonl)
retain every run. These small-field complete-workload checks establish
functional behavior; they do not establish generalization to larger fields.

Primary accounting remains unmeasured:

| Variant | Common operations | S | Ratio to rho | Ratio to floor | Verified targets | Classification |
|---|---|---|---|---|---|---|
| Matched rho reference | null | null | null | null | Unmeasured | Not run in this adapter |
| Preserved pair baseline | null | null | null | null | 72/72 | Reference |
| Opt-in learned selector | null | null | null | null | 72/72 | Engineering capability; no gain established |

The [degree-131 check](degree131-guard.json) initializes the actual public
parameters, checks the model guard and selects pairs for an unseen domain. It
issues zero decomposition queries and builds zero pair entries. Full-size
completion and learned full-size yield remain unmeasured.

## Validation and artifacts

[validation.json](validation.json) records 8 selector tests, 10 existing fixed
workflow tests and 23 Rust CLI/framework/progress tests, all passing. Coverage
includes mixed-signal requirements, planted-control exclusion, orbit leakage,
held-out model selection, bounded fallback after SAT timeout, exact witness
checks, domain guards, invalid models, partial construction, persistence,
completed-result replay and CLI model-path forwarding.

Source hashes are in `baseline-manifest.json` and
`candidate-v2-manifest.json`; the first rejected candidate has its own manifest.
The accepted model and all accepted panels are under `model.json` and
`panels-v2/`. Earlier raw data and source snapshots are retained. See [Publication contents](PUBLICATION.md)
for the scope of the portable evidence and the historical local manifest.
The current runtime implementation matches the frozen accepted candidate.
See `runtime.json` for environment provenance and `manifest.json` for artifact
hashes. The canonical scoreboard adds this result without replacing earlier
measurements.

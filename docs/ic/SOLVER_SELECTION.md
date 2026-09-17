# Learned solver and budget selection

## Restore the research evidence

Large text artifacts are stored as readable chunks so each publication request
fits the automatic review limit. Restore their exact bytes before training on
the saved panels or rerunning an audit:

```sh
python3 research/ic_solver_selection_20260915/restore_evidence.py
```

The script verifies every artifact's SHA-256 hash and refuses to replace a
changed file. Use `--verify` to validate chunks without writing files.

`ic fixed --solver learned --selector-model MODEL.json` uses a small,
cost-sensitive decision tree to choose a decomposition policy for each query.
The default remains `--solver pairs`.

The initial experiment found **no routing headroom**: on all 288 natural
queries, the prepared pair-table policy had the lowest measured cost and the
same verified rank gain as both SAT policies. The fitted model therefore
chooses pairs. This is a working training and selection pipeline, with an
opt-in model; it establishes no algorithmic speedup.

## Use the fitted model

From the repository root:

```sh
cargo build --bin ic
target/debug/ic fixed --params docs/ic/params/k0n9-fixed.json \
  --dir runs/k0n9-learned --solver learned \
  --selector-model research/ic_solver_selection_20260915/model.json \
  --attempts 256 --pair-budget 10000 --query-seconds 0.25 --json
```

The model was trained on degree-9 and degree-23 toy instances. A fingerprint
binds each supported domain to the field, subgroup order, complete signed
factor base, summand count, and query budget. The generator and requested
targets can change without changing that decomposition domain. An unseen
domain or feature range uses pairs and records the reason. In particular,
**ECC2K-130 uses this fallback**; this model supplies no full-size yield
prediction.

Selection happens after bounded pair-table construction. `--pair-budget`
still caps additional construction and `pair_budget` still means incomplete.
The model does not bypass that budget. Resume and relation reuse work as in
the [fixed-parameter workflow](FIXED_PARAMETERS.md).

## Collect labels and train

The panel collector requires the SAT dependencies used by the existing engine
(`python-sat` with CryptoMiniSat and `pycryptosat`). Choose an interpreter that
has them installed:

```sh
/home/ubuntu/crypto-venv/bin/python ecc2k130/codegen/indexcalc_selector_bench.py \
  --params docs/ic/params/k0n9-fixed.json --out runs/selector-panel \
  --queries 128 --seed 1701 --controls 4 --pair-budget 100000 \
  --query-seconds 0.25

/home/ubuntu/crypto-venv/bin/python ecc2k130/codegen/indexcalc_selector.py \
  --data research/ic_solver_selection_20260915/panels-v2/*/queries.jsonl \
  --out runs/selector-model.json
```

The first command collects a single panel. Training usually needs several
panels with a useful mix of natural outcomes; the second command demonstrates
training on the audited mixed panels shipped with this experiment. Output
paths must be new, preserving previous datasets and locked models.

Every query record contains all three policy outcomes on the same target and
the same prior matrix:

| Policy | First call | Fallback when no verified witness returns |
|---|---|---|
| `pairs` | Pair table, full query budget | None |
| `sat-short` | SAT, 10% of query budget | Pair table, remaining 90% |
| `sat-medium` | SAT, 50% of query budget | Pair table, remaining 50% |

SAT's budget limits its native solve. Its existing separate one-second
watchdog bounds the child query; circuit setup, process management and witness
verification add overhead. The sum of fractions is not a hard wall-clock
deadline for the complete policy. Full call times and every timeout remain in
the record. Timeout means unresolved. A final `unsat` comes from completed
exhaustive pair search. Solver errors invalidate a training panel.

Pair and SAT circuit setup are shared and measured once for each panel, so
query costs describe prepared artifacts. Cold complete-workload comparisons
charge construction, model loading, inference, unsuccessful queries, matrix
work, verification, commits, reporting and process exit separately.

## What the model learns

Features available before solving are the query's normal-basis x-coordinate
weight, the weight of its difference from its Frobenius image, current matrix
rank fraction, and whether the query contains the individual target.
No scalar, witness, final solver statistic or planted label is an input.

For each policy, the label is the independently checked matrix rank increase,
and the cost is elapsed nanoseconds including feature extraction and witness
checking. The tree minimizes

```
cost + penalty * (1 - rank_gain)
```

The penalty is the best constant policy's training cost per rank gain. It
prices an unsuccessful query using the observed cost of obtaining useful
progress. This is a cost-sensitive surrogate for useful relations per second;
it is not a proof that the fitted policy minimizes end-to-end work.

The tree has depth at most three and requires at least eight distinct query
orbits per child. Entire signed Frobenius orbits share one train/validation/
confirmation split across seeds, generators and factor bases. Matrix histories
are also kept separate between splits. Planted controls are retained but
excluded from training and natural-yield estimates. Training requires both
natural gains and non-gains. The tree is used only if its held-out validation
loss improves by at least 5%; otherwise the model retains a constant policy.
Confirmation labels do not select the model. Offline selector scores replay
recorded policy costs and exclude inference overhead; the separate complete
workload comparison measures inference overhead in the actual process.

Models use validated JSON trees. Each inference records its features, selected
policy, fallback reason, cost, and model hash in the persistent attempt log.
Every returned relation is checked in the group and its signed-orbit
coefficients are checked again before matrix insertion.

This version implements per-query routing over an existing base. Factor-base
search, a Gröbner adapter, longer decompositions and learned internal SAT
heuristics remain separate experiments.

See the [research report](../../research/ic_solver_selection_20260915/README.md)
for the frozen sources, independent audit, rejected first label dataset,
matched complete-workload results, accounting scope and limitations.

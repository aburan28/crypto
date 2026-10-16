# Bounded natural-target worker for cross-stream N83 pilot bases: source gate

The mixed-orbit example now accepts `run-cross-pilot` for two distinct policy/seed streams from the exact primary-arm 54-row pilot manifest. It replays both compressed archives, reconstructs every signed-Frobenius orbit, forms the exact pool-minus-base residual, and makes **one** bounded natural public-target SAT call with the same frozen first-trial coefficients used by the two retained nested-pair trials. A returned partial is checked against the selected base and residual point domains and independently re-summed to the source equation before it can be reported. A conflict or preflight cap remains `UNKNOWN`; a complete UNSAT trial is recorded separately from a useful partial. This mode has no large-prime graph, rank, linear-algebra, individual-log, or complete-runtime result.

The command shape is:

```text
run-cross-pilot BASE_HASH BASE_MANIFEST BASE_OBJECT POOL_HASH POOL_MANIFEST POOL_OBJECT CORPUS BUDGET_SECONDS CONFLICT_BUDGET MAX_MODELS BASE_SUMMANDS RESIDUAL_SUMMANDS OUTPUT_DIR
```

The worker enforces the same 4-GiB, zero-swap Linux memory cgroup as the nested run modes, records its embedded source digest in `config.json`, checks the public corpus digest, writes an immutable single-call `result.json`, and has an internal wall watchdog. It does not fetch S3 objects or validate a user budget grant. An external, separately authorized supervisor must provide exact archive bytes, a one-CPU/no-network container, charge all launched stages, bind the worker build to reviewed source, and independently audit the result before any retained experiment. Direct invocation is not an approved measurement protocol.

The existing [paired capacity formula](../paired-capacity-screen-20261013/screen.py) gives a useful bound for the [648 primary pilot pair addresses](../cross-stream-pilot-import-source-20261013/census.json). If the exact residual is nonempty, it has at most `K_pool` signed orbits. Evaluating all 648 ordered pairs at that **worst-case** residual size gives:

| Summands | Residual summands | Pairs under both hard caps | Maximum variables | Maximum domain clauses |
| ---: | ---: | ---: | ---: | ---: |
| 5 | 1 | 648 | 97,080 | 269,830 |
| 5 | 2 | 648 | 97,080 | 269,830 |
| 6 | 1 | 648 | 120,696 | 323,796 |
| 6 | 2 | 648 | 120,696 | 323,796 |

Both maxima are below the frozen 150,000-variable and 12,000,000-domain-clause limits. The formula was previously checked against all three counters in the K=600 and K=1,182 retained trial events; this cross-stream calculation is conditional source arithmetic, not an observed SAT allocation or memory measurement. Exact overlap, including whether a pair has an empty residual, still requires archive import.

The input receipts and formula reproduce with `screen.py --check` and `census.py --check` in the linked directories. To reproduce the table, iterate ordered primary objects whose policy or seed differs, call `screen.requirements(K_base, K_base + K_pool, summands, residual_summands)`, and take the maximum `total_variables` and `domain_clauses` for each of the four shapes. `K_base + K_pool` is only a way to pass the worst-case residual size into the nested-pair formula; it does not assert that the two cross-stream bases are nested.

Validation at this source: `cargo test --release --example koblitz_n83_mixed_orbit_worker -- --test-threads=1` passed six tests, including fail-closed same-stream and invalid-budget preflight before archive I/O. The N83 study Python suite passed 63 tests; `rustfmt --check` passed. The full release library suite passed 4,001 tests with 113 ignored on the unchanged library source at the preceding commit. No cross-stream archive was downloaded or imported, no solver call was launched, and no S3 object was written here. The conservative local retained-N83 ledger remains **17,400 seconds**. A minimum-total-runtime base is still unselected.

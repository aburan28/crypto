# n83 F6 wide coordinate screening: exact and faster on planted blocks

The [registered protocol](PROTOCOL.md), [search-status amendment](AMENDMENT_1.md),
and [complete-block comparison amendment](AMENDMENT_2.md) are the rules
for these results. This is an **exploratory three-summand F6 algebraic
block**, not a complete n83 F6-IC candidate or an ordinary relation.
The pinned K0 curve is `icv1-f2m83-tm6151469093347-debefd74`; the
dimension-12 source subspace has 4,057 curve points and 2,029 distinct
valid abscissae. Cofactor-four projection has 4,054 usable subgroup
points. All probes use planted source indices `[0,2,4]`.

## Change and correctness

`valid_coordinates` previously solved the curve equation at every
subspace abscissa on every call. It now verifies that `index_of` indexes
the supplied factor base, takes allowed abscissae from its existing
points, and builds the sorted coordinate list with field-bit XORs. A
mismatched index retains the old curve-lift path. An exact test compared
the new and old lists at n=9, 17 and 83, with restricted indices,
a changed basis and a malformed index. The n83 list contained exactly
2,029 coordinates.

The first release suite [failed](wide_tests_initial_failure.log.gz):
the new test used an unavailable n17 fixture, and an existing suffix
test hit an intermittent incomplete `Infinity` remainder. After fixing
the fixture, the indexed-list [reference test](indexed_reference_test.log.gz)
passed, while one of ten direct suffix replays still failed. The search
had treated an unsupported branch as a global stop. Its internal status
now lets other branches continue; an unsupported **no-witness** result
still sets public `exhausted=true`, so callers cannot call it a proved
miss. A verified witness clears both incomplete flags. The final
[nine release wide tests](wide_tests_final.log.gz) passed, including
the legacy unsupported-status test. The suffix test passed
[20 of 20 direct replays](suffix_final_20.status.tsv); the combined
[raw log](suffix_final_20.log.gz) and initial ten replay logs remain
in this directory. All failed runs were preserved.

## Paired local measurements

Five alternating pairs used frozen release executables, the same exact
source target and a 120-second cap per process. All 20 process exits
were zero and every stderr file was empty. Setup is measured separately
from the query interval; process launch is outside both. The host was
not isolated, so these are exploratory stage diagnostics and observed
paired ranges, not controlled CPU speedup claims.

| Frozen workload and outcome | Baseline median query | Candidate median query | Query ratio | Setup + query median, baseline → candidate | Peak RSS median, baseline → candidate |
| --- | ---: | ---: | ---: | ---: | ---: |
| One-node root; both `exhausted`, no witness | 66.228 ms | 17.878 ms | 3.704× | 117.075 → 68.768 ms (1.702×) | 21.053 → 19.923 MB |
| Complete planted block; both found `[0,2,4]` | 72.437 ms | 18.754 ms | 3.862× | 125.274 → 74.557 ms (1.680×) | 20.562 → 19.530 MB |

Every root run had two nodes, two reductions and maximum matrix
166 × 15,259. Every complete planted run had four nodes, five
reductions, one leaf, the same matrix dimensions, and exact source-sum
and cofactor-four projection replay. Root paired query ratios ranged
2.635–3.985× (paired median 3.732×); complete planted paired query
ratios ranged 2.694–4.306× (paired median 4.018×). All five paired
planted ratios exceeded the registered 2× local gate. Candidate median
RSS was lower on both workloads, within the registered limit. The
individual [root runs](status.tsv) and [complete-block runs](witness_status.tsv)
are retained as JSON lines beside their status and stderr files. The
single candidate [pre-pair witness](witness_probe.jsonl) also found and
replayed `[0,2,4]` in four nodes.

The setup-plus-query figures are a **cold probe proxy**, not cold IC
time: they include curve, source-base, field-table and index-map setup
but omit relation collection, matrix solving, target descent and scalar
recovery. The query figures are a planted block slice, not a natural
single-target online IC interval. Neither table row establishes a 2×
complete F6, F4/F5 crossover, end-to-end IC speedup, or rho crossover.
The incomplete pipeline has no `IC1` candidate ID or verified run ID.

## Reproduction and custody

The final tests ran `cargo test --offline --release --lib
cryptanalysis::wide_groebner::tests::`; the [probe build log](probe_build.log.gz)
covers `cargo build --offline --release --example
f6_n83_wide_block_probe --example f6_n83_wide_lazy_witness_probe`.
Cargo used `CARGO_TARGET_DIR=/Volumes/SSD990/crypto/worktrees/f6-n83-geometric-20261004/target`
and `TMPDIR=/Volumes/SSD990/crypto/worktrees/f6-n83-batch-size-20261004/target/tmp`.
Run `sh research/f6_n83_valid_coordinates_20261005/run.sh` and
`sh research/f6_n83_valid_coordinates_20261005/run_witness.sh` from the
repository root against the frozen binaries in `target/bench_bins`.
The [baseline witness rebuild log](witness_baseline_rebuild.log.gz)
records rebuilding the exact parent source in this same worktree before
the complete-block pairs.

Source commit `b25e9f4cd`; baseline parent
`128cf786a60da4f4289d8680092916991cf90f8a`. SHA-256:

| Artifact | Digest |
| --- | --- |
| Baseline/candidate root binaries | `16bfcd5408a8c736269b9a2d1359375f8bd3eb2880c194bf11a0333dfd7db22f` / `a6e2ae0bb26dfaac7c2653196bdec683c05796493690336ba5a34996028e9c12` |
| Baseline/candidate complete-block binaries | `8d79eebbb93634a7a4faa3f58ef3cb1827797d079abe22ef35aedb1dd3cc7cd2` / `55a25c8ccf25e5c9184990160a586ebe41ec3e1f5b835f5dd5715ca097f23eeb` |
| Solver source | `af247b414dd315ab5371322226d88a6717d95b7dffbf2152b55edb09fdbc774f` |
| Curve/factor-base source | `1c98b492fafe80ae3a9dbea5400e7fceb114dec7b9a548e4d98f33a691f05b5d` |
| Root/witness probe sources | `fde4a969957de325ab0803adba429382165b119dd97b70ecc0c0d0763166659c` / `109c24b6cfcf42b20b016190df04e743299409ea7fa21ed0bc7832551ad78e66` |
| Protocol, Amendments 1 and 2 | `d9ba2ee203324e6c726362eec81a415abb9982c0cb091805e93f6466868e61d7` / `99a25aa504aa5d50547b48dce22c503769840cfaf1d53c2a4a0aadb30f49ba07` / `0fdf058f91e2415a5953270f0625d9b742ded6ba660ce190a971b4aaca44fccf` |

Physical Apple M4 Pro, 14 cores (10 performance, four efficiency),
arm64 macOS 26.6, Rust 1.93.1 release builds. No host-level CPU
isolation receipt exists; the memory numbers are process `getrusage`
peaks. There is no x86 or GPU performance claim.

The remaining n83 method gate is an ordinary verified relation from
eight or more source summands, charging all failed attempts. The
128-variable block solves three summands; it does not solve the full
high-arity equation. A completed relation matrix, fresh target DLP,
same-point one-target rho reference and isolated-host timing are still
needed before any end-to-end speedup can be stated.

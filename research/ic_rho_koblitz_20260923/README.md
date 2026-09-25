# The Koblitz references: batch rho, the step and the stored pair, priced

Evidence for §19 of
[`research/notes/index-calculus/RESEARCH_IC_BOUNDARY_LEDGER.md`](../notes/index-calculus/RESEARCH_IC_BOUNDARY_LEDGER.md).
The Koblitz collection thread quoted every `vs rho` as the cost of 32 targets
solved together, per target, divided by rho solving one. It also counted each
rho step and each stored pair of its folded pair table as one addition. The
repository's multi-target rule (`src/ecc_safety.rs::check_multi_target_margin`,
Kuhn–Struik) makes batch rho at the same `k` the reference for such a figure.
§19.1 declared the boundaries, three measurements, three targets and the
abandon conditions, and was committed before anything here ran.

## What was measured

| part | what | code |
|:--|:--|:--|
| (a) batch rho | `k` targets solved in sequence by the tuned walk of §18 on the signed-Frobenius classes (`A = 2n`), jumps in `G` only, one distinguished-point table per batch, so a later target can finish on an earlier one's trail; fresh targets every batch, every logarithm checked against the planted one | `ic_boundary::rho_batch_with`, `signed_frobenius_rho_batch`; `ic rho --batch-koblitz` |
| (b) the step | in the thread's unit, one batched addition (`add_many` over 1,024 points): a canonical step (the unit plus `SignedFrobeniusClasses::canon`) and Bailey et al.'s step (the unit plus the normal-basis weight and two table-applied Frobenius powers) | `examples/koblitz_reference_prices.rs` |
| (c) the stored pair | the folded table's whole build over the pairs it stores, in the same unit, on one thread | the same |

## Layout

| path | what | ledger |
|:--|:--|:--|
| `batch/k0n41.json`, `k0n53.json` | `K_0/GF(2^41)` and `K_0/GF(2^53)`, `k = 1, 4, 16, 32`, 16 batches each | §19.3 |
| `batch/k0n61.json` | `K_0/GF(2^61)` (the `n = 61` panel's curve), `k = 1, 32`, 8 batches | §19.3 |
| `batch/*.stderr.log` | each batch's progress line and the summary table | |
| `prices/n41-F15744.json` | (b) and (c) on the thread's `n = 41` base, `\|F\| = 15,744`: the declared price | §19.4 |
| `prices/n53.json`, `n61.json` | (b) on the other two curves | §19.4 |
| `prices/n41-F5248.json`, `n41-F16400.json`, `n53-F15264.json`, `n61-F*.json` | (c) on the base of every other quoted figure; not declared, so that no base carries another's price. `n41-F16400` was added after the first pass and before any analysis | §19.4 |
| `analyse.py`, `analysis.json` | grades the three targets and re-reads all 37 quoted figures from their frozen `docs/ic/runs` files | §19.5 |
| `provenance.txt` | the commit (clean tree), the binaries' hashes, the host, the times | |
| `run.sh` | every command, in order | |

Each price file also carries two diagnostics that §19.1 did not declare. One
is `add_pairwise`, the batched addition that 1,024 independent walks actually
run. The other is a lone affine addition, the unit of the three-regime ledger's
single-target rows.

## Results, in brief

- **All three declared targets met.** All 1,960 targets were recovered and
  verified. Per-target cost at `k = 32` over `k = 1` is `0.210`, `0.239` and
  `0.190`, inside the declared `[0.09, 0.35]` and each interval covers the
  Kuhn–Struik `0.199`. All 37 figures were re-read.
- **Batch rho at `k = 32`** costs `0.0300`, `0.0252` and `0.0219` a target.
- **A canonical step costs `2.74–2.92` units, Bailey's `1.38–1.41`, a stored
  pair `4.5–6.1`** (`6.14` on the declared base). The folded build makes two
  passes over every row.
- **The closest figure, `1.17×`, reads `6.51×`** against batch rho as counted,
  and `6.38×` with rho's step and the build both priced. It reads `12.6×` if
  Bailey's step is taken with the canonical walk's count, which is a model.
  One target alone costs `43.7×` rho.
- **Class: accounting.** No index-calculus count changed.

## Reproducing

```
cargo build --release --bin ic --example koblitz_reference_prices
sh research/ic_rho_koblitz_20260923/run.sh
python3 research/ic_rho_koblitz_20260923/analyse.py
```

`ic` never overwrites a report, so move the frozen files aside first. The
price measurements are timings: run them alone on the host, as `run.sh` does,
before the batch runs start.

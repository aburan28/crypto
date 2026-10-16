# Amendment 4 (2026-10-15, before any data for these arms)

## Why

Amendments 1–2 showed CADO-NFS's filter leaves a 5–8× smaller matrix than
ours. Reading its logs located the gap in `purge` (whole-clique removal over
up to 50 steps), not only in `merge`. This change ports both into
`koblitz_sparse_la::filter_relations` as `FilterStrategy::Cado` (purge: Ω₃₁
clique removal; merge: Markowitz order, minimum-spanning-tree combination,
density stop) and ships `FilterOptions::cado()` (merge columns ≤ 32, no row
cap, target density 100), now the default of `ic logs`. The legacy filter
stays the library default so frozen job configurations replay unchanged.

**Disclosure.** The shipped settings were chosen from development probes on
seed 1 at D = 65,536 and 262,144 (`+1` coefficients; target density 0
against 100). Those probes are not data. The density sweep below is in the
protocol to show the sensitivity of that choice.

## Arms (m = 3, seeds 1–3, four threads, 1,800 s timeout, interleaved by
arm within each (D, seed), one fresh session so every arm sees the same
host state)

| arm | coefficients | filter | Krylov |
|:--|:--|:--|:--|
| `legacy-bw` | ±1 | `FilterOptions::default()` | block Wiedemann |
| `cado-bw` | ±1 | `FilterOptions::cado()` | block Wiedemann |
| `cado-bl` | ±1 | `FilterOptions::cado()` | block Lanczos |
| `legacy-bw-pos` | +1 | default | block Wiedemann |
| `cado-bw-pos` | +1 | `cado()` | block Wiedemann |
| `cado-bw-pos-d0` / `-d50` / `-d170` | +1 | `cado()` with target density 0 / 50 / 170 | block Wiedemann |
| `cado-nfs` | +1 | CADO-NFS purge + merge-dl + replay-dl + bwc (`run_cado_full_fb.sh`) | CADO bwc |

D = 1,024 … 262,144 in powers of two. Every arm verifies all D logs
against the planted vector. Headline ratios are per (D, seed) within this
session: `legacy-bw / cado-bw` and `cado-nfs / cado-bw-pos`.

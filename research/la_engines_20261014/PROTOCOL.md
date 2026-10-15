# Protocol: index-calculus linear-algebra engines, ours against CADO-NFS `bwc`

Frozen 2026-10-14, before any measurement in this study. Changes after data
are listed under "Deviations" in the report, never edited in here.

## Questions

1. **Engine parity.** On an identical square relation core over GF(p), how
   does our block Wiedemann (`koblitz_sparse_la`) compare in wall time with
   CADO-NFS `bwc` (block Wiedemann, GF(p), inhomogeneous right nullspace)?
2. **Our own engines.** Across factor-base sizes D, how do the four solvers
   we ship compare: `IncrementalGauss` (ecbench default), `StructuredGauss`,
   filter + block Wiedemann, filter + block Lanczos?
3. **End to end.** In a real `ic logs` precompute (relation collection,
   then the linear algebra, then the certification), what fraction of the
   wall time is linear algebra, for each `--linear-algebra` mode?

## Instrument

- `examples/la_engine_bench.rs` (Rust). Synthetic systems: D columns, rows
  of m = `--weight` summands at uniform columns with coefficients ±1
  (generic prime-curve shape), drawn until every column occurs, then 32
  more rows. A planted uniform logarithm vector over
  p = 2^62 − 57 fixes each right-hand side. Deterministic per
  (D, m, seed).
- The CADO core: singleton + clique filtering with merge disabled,
  trimmed to square (`unmerged_core`), so coefficients stay small integers
  that the `bwc` binary format can carry. `core-bw`/`core-bl` and
  `cado-bwc` see exactly this matrix and right-hand side.
- CADO-NFS commit `692ecb7e62f0f3bdab88ee44cc60b8ded0ea1a1b` (2026-10-02),
  Release build, GF(p) layers `p1;p3`, `bwc.pl :complete` with
  `prime=p nrhs=1 nullspace=right mpi=1x1`, timed from `bwc.pl` start to
  exit (includes its balancing, `prep`, `secure`, `krylov`, `lingen`,
  `mksol`, `gather`; excludes our export and `mf_scan2`).
- Every arm's answer is checked natively against the planted vector
  (`verified`). An unverified run is reported, never dropped.

## Arms and grid

| arm | what | D |
|:--|:--|:--|
| `igauss` | dense incremental RREF, rows fed until rank D | 1k–8k (stop when one run > 600 s) |
| `sgauss` | sparse incremental RREF | 1k–16k (same stop rule) |
| `full-bw` | filter + block Wiedemann 4×4 + reconstruct | 1k–256k |
| `full-bl` | filter + block Lanczos 4 + reconstruct | 1k–256k |
| `core-bw` | ours, block Wiedemann 4×4 on the unmerged core | 1k–256k |
| `cado-bwc` | CADO `bwc` m = n = 4 on the same core | 1k–256k |

D runs in powers of two from 1,024. m = 3 throughout; m = 2 and m = 4 at
D = 16,384 only. Seeds 1, 2, 3 per cell. Threads: 4 (Rayon pool / CADO
`thr=2x2`), and a 1-thread pass (`RAYON_NUM_THREADS=1` / `thr=1x1`) for
`core-bw` and `cado-bwc` at D ≤ 65,536. Per-run timeout 1,800 s; a timeout
ends that arm's sweep upward and is reported.

Order: cells interleaved by arm within each (D, seed) to spread drift.

## Controls

- **A/A:** `core-bw` repeated on seed 1 at every D (same matrix twice);
  the spread bounds what a single-run difference can mean.
- **Positive control:** every arm recovers the same planted vector, so a
  "fast" arm that solves nothing cannot pass.

## End to end (question 3)

`target/release/ic logs` on Koblitz curves at the largest configurations
that complete within 1,800 s on this host, with `--linear-algebra sparse`
(block Wiedemann and block Lanczos) and `--linear-algebra dense`, reading
`collection_seconds` and `linear_algebra_seconds` from the report. Three
seeds per configuration.

## Reporting

Wall seconds, medians over seeds with min–max, per arm and D; ratios
`cado-bwc / core-bw` per (D, seed). Host: 4 vCPU Intel Xeon @ 2.80 GHz
cloud VM, shared, no core reservation (isolation L0): wall times are
descriptive, not admitted ecbench measurements. No ECDLP speedup claim
follows from any of this; the study prices one phase.

# Joux–Vitse `(n − 1)`-decomposition at `n = 5`: measured end to end against rho

2026-10-15.  Report with figures: [`REPORT.pdf`](REPORT.pdf) (source
[`REPORT.typ`](REPORT.typ), which reads [`summary.json`](summary.json)).  Run
report: [`reports/run-20261015T0518Z/`](reports/run-20261015T0518Z/).  Cost
model and literature: [`MODEL.md`](MODEL.md).  Protocol frozen before data:
[`PROTOCOL.md`](PROTOCOL.md), with [`PROTOCOL_ADDENDUM.md`](PROTOCOL_ADDENDUM.md).
Code: `src/cryptanalysis/joux_vitse.rs`, `examples/joux_vitse_ladder.rs`,
`examples/joux_vitse_summary.rs`.  Open obstruction:
[`problems/k5-joux-vitse-trial-cost-and-fpk-ecbench.md`](../../problems/k5-joux-vitse-trial-cost-and-fpk-ecbench.md).

**Class: a measurement of a published route (reproduction), plus an
accounting correction to the scoreboard's `k = 5` plain constant. No
advance.**

## Verdict

No. The `(n − 1)`-decomposition does not give this repository a route where
measured IC cost approaches rho at sizes it can run.

- **Complete DLPs.** IC is `4.5·10⁷`–`1.8·10⁸` times rho, measured on complete,
  verified DLPs at `q = 11, 31` against rho run on the same target, in the
  same unit (`F_p` multiplications).
- **Trial cost.** One trial costs `7.5·10⁹` `F_p` multiplications on average,
  flat in `q`.
- **Trials and relations.** About `25 q` trials buy one relation, and `q/2`
  relations are needed.
- **Crossover.** On these constants the method meets rho at `q ≈ 2^{58.2}`, a
  group of about `2^{291}` (extrapolation). Vitse's thesis states `2^{60}`.
- **At runnable sizes.** The model gives about `10⁷` times rho at
  `q = 2^{10}` and `2·10⁶` at `q = 2^{16}`.

### Why it fails

The relation search costs `≈ 12 q² · C_trial`, against rho's `≈ 0.89 q^{5/2}`
group operations. The ratio therefore falls only as `q^{−1/2}`, from a
constant set by `C_trial / c_add ≈ 4·10⁷`.

Large primes cannot help. There is no linear algebra to rebalance, and a
smaller factor base raises the trials per relation as its fourth power
(MODEL.md).

## Results

### One trial

There are 12 rungs (`p = 11 … 241`), each with 2 curves × 32 F4 trials,
every one cross-checked. On top of those, `500 p` meet-in-the-middle trials
per curve measure the decomposition rate. Source: `data/ladder/`,
`summary.json`.

| `p` | `C_trial` (mean, `F_p` mult.) | trials / relation (Wilson 95 %) | rate · 24p | `C_rel` |
|--:|--:|--:|--:|--:|
| 11 | 9.21·10⁹ | 503 (332–761) | 0.52 | 4.63·10¹² |
| 31 | 9.17·10⁹ | 1,002 (706–1,422) | 0.74 | 9.19·10¹² |
| 41 | 7.63·10⁹ | 893 (669–1,190) | 1.10 | 6.81·10¹² |
| 61 | 8.54·10⁹ | 939 (737–1,197) | 1.56 | 8.02·10¹² |
| 71 | 6.51·10⁹ | 2,030 (1,460–2,823) | 0.84 | 1.32·10¹³ |
| 101 | 8.72·10⁹ | 2,021 (1,533–2,664) | 1.20 | 1.76·10¹³ |
| 131 | 7.11·10⁹ | 3,121 (2,309–4,218) | 1.01 | 2.22·10¹³ |
| 151 | 6.22·10⁹ | 3,684 (2,716–4,998) | 0.98 | 2.29·10¹³ |
| 181 | 8.91·10⁹ | 3,772 (2,845–5,001) | 1.15 | 3.36·10¹³ |
| 191 | 6.08·10⁹ | 7,349 (5,015–10,768) | 0.62 | 4.47·10¹³ |
| 211 | 5.91·10⁹ | 4,690 (3,506–6,275) | 1.08 | 2.77·10¹³ |
| 241 | 6.33·10⁹ | 6,027 (4,426–8,206) | 0.96 | 3.82·10¹³ |

**Fits** (least squares on logs, 12 rungs):

| quantity | fitted exponent of `p` | predicted |
|:--|:--|:--|
| `C_trial` | `−0.12 ± 0.05` | `0` |
| trials per relation | `0.87 ± 0.09` | `1` |
| `C_rel` | `0.75 ± 0.08` | `1` |

The trials-per-relation slope sits below 1 because `|F| / p` varies from
curve to curve at small `p`. Above `p = 41` the rate matches `1/(24p)`:
rate · 24p = 0.62–1.20.

**Checks failed: 0.** On every F4 trial, the set of decompositions into four
distinct base points equals the independent meet-in-the-middle set, and
every relation re-adds to `R` in the group.

**Replayed trials cost one of three levels:** `≈ 5.7`, `8.2` and
`11.1·10⁹`. The last equals a full F4 run (`1.07–1.10·10¹⁰`). The suspected
cause is that `f4_fp` finishes pairs still pending after the recorded steps
as a full run. This is not verified per trial.

### Complete DLPs

Each IC run is paired with rho on the same `(P, Q)`. Source: `data/dlp/`.

| run | `ℓ` | IC trials | relations | IC online | rho ops | IC / rho | IC / expected walk | answer |
|:--|--:|--:|--:|--:|--:|--:|--:|:--|
| `p11_s1` | 160,637 | 1,912 | 6 | 1.74·10¹³ | 2,030 | 4.5·10⁷ | 1.8·10⁸ | both verified |
| `p11_s2` | 4,493 | 2,920 | 7 | 2.67·10¹³ | 800 | 1.8·10⁸ | 1.7·10⁹ | both verified |
| `p11_s3` | 14,669 | stopped | 0 | — | 1,054 | — | — | rho only |
| `p31_s1` | 28,622,519 | 6,688 | 19 | 5.31·10¹³ | 4,092 | 6.8·10⁷ | 4.2·10⁷ | both verified |

- **IC online** is in `F_p` multiplications. It covers target generation,
  the descended systems, the F4 trials (replays and full fallbacks), solving,
  splitting, re-adding, and linear algebra (multiplications modulo `ℓ`).
- **Set-up is excluded from IC online and reported separately.** This is
  target-independent work: the symmetrised `S₅` (`1.87·10⁹`) and trace
  recording (`≈ 10¹⁰`).
- **Rho's counts include its 20 multipliers' set-up.** That set-up dominates
  at the `p = 11` orders. "IC / expected walk" prices rho as its expected
  walk alone, `√(πℓ/2)` operations at the measured `c_add`.
- **`p11_s3` was stopped.** It drew a curve with a one-point factor base
  (`data/dlp/p11_s3.STOPPED.md`). No other seed was substituted.

Of 672 consistent systems in the three complete runs, 702 of the 740
`F_p`-solutions found did not split. Those are decompositions over an
extension, not relations. This is the expected `≈ 23/24`.

### Attribution

The committed `jv_quintic` plain test, re-run unchanged on today's `f4_fp`
(`data/jv_quintic_rerun/`), costs `C″ = 1.104·10¹⁰`. The scoreboard says
`3.0·10¹⁰`, measured 2026-09-29.

- **The `2.7×` from `3.0·10¹⁰` to `1.1·10¹⁰`** is the engine's own progress
  (Gebauer–Möller criteria, packed monomials). It is not credited here.
- **Trace replay** is credited with the rest: `1.1·10¹⁰ → 7.5·10⁹` mean
  (`5.7·10⁹` at its best level).
- **Plain-route crossover:** on the corrected constant, the scoreboard's
  `2^{302}` becomes `≈ 2^{291}`.

## Scope and limits

- **Survey miss.** The repository already measured the `k = 5` test as a
  stage diagnostic: the plain and 2-torsion Edwards variants (`jv_quintic.rs`,
  `jv_quintic_edwards.rs`, `RESEARCH_K5_TORSION_JOUX_VITSE.md`). This study's
  survey missed it before building, and `joux_vitse.rs` duplicates the plain
  construction. What this study adds:
  - the first complete `k = 5` DLPs;
  - trace replay of the trial;
  - the rate ladder with its oracle;
  - the corrected constant.
- **Better per-trial constant elsewhere.** The 2-torsion Edwards test
  (`7.3·10⁷`) is a better constant than anything here. Its crossover,
  `2^{233}`, is likewise far out of reach.
- **Correction to the premise: the practical `F_{p⁶}` result.** Joux–Vitse's
  practical `F_{p⁶}` result is cover and decomposition (thesis ch. 8), not
  `(n − 1)`-decomposition. The repository measures that route on a weak
  class, below rho from `p ≈ 430`; see the rho-parity ledger.
- **Correction to the premise: the practical `F_{q⁵}` result.** Their
  practical `F_{q⁵}` result is oracle-assisted static Diffie–Hellman, not a
  DLP (MODEL.md).
- **The PR #1863 comparison** (`n = 4` D₄ oracle, `C₄ = 6.9·10⁷`) is
  per-relation on a different group, `q⁴` against `q⁵`. #1863 is open and
  stacked on the unmerged #1792.
- **`n = 6`: not attempted.** The symmetrised `S₆` needs about 20 GB against
  this host's 15 GB, about `1.3·10¹⁴` multiplications, and a trial with
  Macaulay bound 91. This is derived, a resource limit.
- **No ecbench run.** `main` has no `F_{p^k}` instance; `Instance::Fp4` is in
  #1792. The problem record holds what would unblock it.
- **Curves are not registered.** Their ICV1 identities are computed in
  [`curves.json`](curves.json), but the registry builder harvests no
  extension-field JSON.
- **Host.** Claude Code cloud container, 4 vCPU Xeon @ 2.10 GHz, 15 GB,
  isolation none. Counts are the metric; wall time is a practicality note.
  - The frozen runs used commit `1e8592da0`. After the runs, two clippy-only edits replaced
    `n % d == 0` by `n.is_multiple_of(d)` (same behaviour).
  - A container restart killed the first launch. Its partial records are in
    `data_killed_20261015T0518Z/`, and the restarted runs reproduced their
    relations bit for bit.

## Reproduce

```sh
cargo build --release --example joux_vitse_ladder --example joux_vitse_summary
sh research/joux_vitse_n1_20261015/run.sh ladder   # and dlp11, dlp31
./target/release/examples/joux_vitse_summary research/joux_vitse_n1_20261015/data \
    > research/joux_vitse_n1_20261015/summary.json
python3 research/joux_vitse_n1_20261015/figures/plot.py   # layout only
```

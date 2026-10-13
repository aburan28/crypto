# Double large prime on the Koblitz pair-table pipeline (2026-10-12)

**Question.** `docs/ic/LEAD_LADDER_METHODOLOGY_20261010.md` §7 item 4 closed the
single-large-prime variation on the Frobenius-folded table by analysis: with the
residual collector as the `|D| = r` limit, the break-even is
`|F| < 2.4·(K r)^{1/4}`, and "a bounded `D` lowers both the entry rate and the
collision rate by the same factor". That analysis was for **one** large prime
matched by *membership* of the residual in `D`. This directory wires one- and
two-large-prime partial relations into the production path and measures whether
the double-large-prime graph — whose partials are found by a *pair table over
`D`*, not by membership — changes the conclusion.

Public synthetic instances (`K_0: y² + xy = x³ + 1` over `F_{2^n}`, planted
logarithms), fail-closed evidence (every partial and every merged row replayed as
a group equation twice: in the exact eliminator and in the counted fast group),
operation counts as the metric. No key-recovery, sub-rho or deployed-curve claim
is made or implied; everything below is a constant on a toy rung inside the
pair-table family.

## What was wired (the seam)

* `ic_boundary::Relation` and `collect_and_solve_relations`: the shared relation
  loop now accepts, besides summand indices, rows with explicit column
  coefficients and their own `(a, b)` (`Relation::Row`). The classic entry point
  `collect_and_solve_with_completion` wraps it and produces bit-identical
  counters (tested).
* `DecompositionOracle::decompose_relations` (default: `decompose` wrapped) and
  `relation_counters` (default: none): the framework runner calls these; every
  existing oracle is unchanged.
* `ic_framework::large_prime::LargePrimeCollector`: the three-summand
  `mitm-frobenius` loop with a large-prime domain `D` = the next
  `lp_domain_factor · |F|` signed points of the compact-orbit raw-x scan after the
  base, a second Frobenius pair table over `F ∪ D` without the `F × F` pairs (and
  without `D × D` for `large_primes = 1`), partial shapes `R = P_i + L`
  (orbit-key lookup), `R = P_i + P_j + L`, `R = P_i + L_1 + L_2`, online merging
  through `koblitz_large_prime::WideLargePrimeEliminator` (exact `±λ^t`
  signed-orbit coefficients modulo `r`), and an independent replay of each
  merged row in the counted group before it reaches the matrix. One complete
  relation per trial, as the classic oracle.
* `FrobeniusPairTable::build_filtered` / `orbit_key` and
  `WideLargePrimeEliminator::new_owned` / `target`: the two small API additions
  the collector needed.
* `ic bench --large-primes {0,1,2} --lp-domain-factor c`, or per configuration
  `mitm-frobenius:m=3,large_primes=2,lp_domain_factor=10`.

Tests: `cargo test --release --lib ic_framework::large_prime` (planted one- and
two-large-prime merges with exact expected rows, wrong coefficient rejected,
`large_primes = 0` identical to the classic loop, end-to-end recovery with 1 and
2 large primes) and `cargo test --release --lib koblitz_large_prime`.

## Reproduce

```sh
export CARGO_TARGET_DIR=$PWD/target
cargo build --release --bin ic
cargo test --release --lib ic_framework::large_prime      # 3 tests
cargo test --release --lib koblitz_large_prime            # 4 tests
cd research/ic_techniques_20261012/double_large_prime
python3 make_sweeps.py 23 31 41                           # sweeps/dlp_n*_c*.json + the n=41 control
./run.sh 23                                               # 3 sweeps, ~3 min wall here
./run.sh 31                                               # 3 sweeps, ~3 min wall here
./run.sh 41                                               # c10 took ~50 min wall on a load-100 host; c5 and c20 not run (below)
SWEEP_PATTERN=control_envelope ./run.sh 41                # the equal-envelope control, ~5 min wall
python3 summarise.py runs/20261012T2105Z_n23 runs/20261012T2305Z_n31 runs/20261012T2315Z_n41 runs/20261013T0030Z_n41_control
```

Each run directory holds the exact command (`cmd_*.txt`), stdout, the commit
(`39b6d66191`, the branch's `origin/main` base; the worktree's uncommitted
changes are the ones this commit lands), the binary's SHA-256
(`1aceeb36…a0605e`, rebuilt from the final source and found identical) and the
host load at start and end. `results.json` is the cell table the summariser
wrote from the run directories listed in its `sources`; the table below is its
markdown rendering. Run-directory stamps are launch labels (the control's
`20261013T0030Z` label was written by hand and is ahead of its actual launch,
2026-10-12 23:37 UTC, which the `host_load_at_start.txt` local time shows).

Instances: `ic bench --koblitz-a 0` at degree `n` is `K_0: y² + xy = x³ + 1`
over `F_{2^n}` with the prime-order subgroup of the registry curve
`icv1-f2m23-t5197-69e76b73` (`r = 2,095,853`, 21.0 bits),
`icv1-f2m31-tm90707-c95f16f5` (`r = 1,439,393`, **20.5 bits** — the group
order is `4 · 373 · r`, so the n = 31 rung is not larger than n = 23) and
`icv1-f2m41-tm2308219-7f48b14a` (`r = 549,756,390,943`, 39.0 bits, cofactor
4). Targets: `walk`, seed `20261012`, 8 planted targets per cell at 5
columns and 4 at 10 and 20; one target pins after about `columns + 1`
complete relations, so a cell holds 28–72 complete relations (≥ 32 was
asked; the 10- and 20-column cells at 4 targets hold 28–44 and 53–72). The
host was a shared 14-core Apple-silicon machine at load average 67–155
throughout; wall time is a practicality note only.

## Results

Columns: `trials / rel.` is total trials over total complete relations of the
cell (direct hits and merged rows together); `partials fed` are the partial
relations the collector handed the exact eliminator, `stored at end` the
pivots still open when the last target pinned; `merge yield` is merged
complete rows per partial fed and `merged share` the fraction of complete
relations that came from merges rather than direct hits; `rel. adds / rel.`
is the relation phase's additions and doublings per complete relation
(probe additions, walk steps and the replay scalar multiplications of merged
rows); `setup adds / target` is the oracle's table build (classic table, and
the `F ∪ D` table where present); `LP table entries` is the `F ∪ D` table's
size; `row ops / target` the incremental Gauss work; `total GAE / target` the
cold one-target total (base build + oracle setup + relation phase + linear
algebra + recovery check) in group-addition equivalents; `S` is total over
`√r`, `S / rho` the ratio to the counted signed-Frobenius rho reference run
in the same report (8 runs per sweep; `S_rho` = 1.08 at n = 23, 0.91 at
n = 31, 0.168 at n = 41); `rejected / capped / replay fail` are the
eliminator's fail-closed counters (a partial whose point equation does not
hold, a merge that hit a limit, a merged row that failed the counted-group
replay). `LP = 0` is the unmodified `mitm-frobenius:m=3` oracle; the `cols =
40` and `110` rows at n = 41 are the equal-envelope control (the classic
oracle on a base as large as `F ∪ D` of the `|D|/|F| = 3` and `10` cells).

| n | cols | \|F\| | LP | \|D\|/\|F\| | targets | trials | complete rel. | trials / rel. | partials fed | stored at end | merge yield | merged share | rel. adds / rel. | setup adds / target | LP table entries | LA rows / target | row ops / target | total GAE / target | S | S / rho | rejected / capped / replay fail | verified |
|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|--:|:--|
| 23 | 5 | 230 | 0 | — | 8 | 78 | 48 | 1.62 | 0 | 0 | — | 0 | 361 | 690 | 0 | 6 | 80.2 | 3242 | 2.24 | 2.07 | 0 / 0 / 0 | yes |
| 23 | 5 | 230 | 1 | 3 | 8 | 46 | 46 | 1 | 112 | 66 | 0.295 | 0.717 | 304 | 4460 | 3311 | 5.75 | 86.4 | 7173 | 4.95 | 4.57 | 0 / 0 / 0 | yes |
| 23 | 5 | 230 | 1 | 10 | 8 | 46 | 46 | 1 | 221 | 157 | 0.167 | 0.804 | 305 | 13070 | 10188 | 5.75 | 86.1 | 16567 | 11.4 | 10.6 | 0 / 0 / 0 | yes |
| 23 | 5 | 230 | 1 | 30 | 8 | 46 | 46 | 1 | 443 | 348 | 0.0903 | 0.87 | 313 | 37670 | 24164 | 5.75 | 83.8 | 43438 | 30 | 27.7 | 0 / 0 / 0 | yes |
| 23 | 5 | 230 | 2 | 3 | 8 | 40 | 40 | 1 | 130 | 89 | 0.277 | 0.9 | 331 | 9980 | 7796 | 5 | 66.1 | 12592 | 8.7 | 8.02 | 0 / 0 / 0 | yes |
| 23 | 5 | 230 | 2 | 10 | 8 | 42 | 42 | 1 | 295 | 251 | 0.125 | 0.881 | 316 | 71720 | 35545 | 5.25 | 74.1 | 75114 | 51.9 | 47.9 | 0 / 0 / 0 | yes |
| 23 | 5 | 230 | 2 | 30 | 8 | 43 | 43 | 1 | 445 | 388 | 0.0854 | 0.884 | 315 | 5.59e+05 | 45562 | 5.38 | 74.6 | 5.64e+05 | 390 | 360 | 0 / 0 / 0 | yes |
| 23 | 10 | 460 | 0 | — | 4 | 37 | 37 | 1 | 0 | 0 | — | 0 | 127 | 2530 | 0 | 9.25 | 185 | 4381 | 3.03 | 2.89 | 0 / 0 / 0 | yes |
| 23 | 10 | 460 | 1 | 3 | 4 | 39 | 39 | 1 | 94 | 70 | 0.223 | 0.538 | 187 | 16970 | 11908 | 9.75 | 176 | 20701 | 14.3 | 13.7 | 0 / 0 / 0 | yes |
| 23 | 10 | 460 | 1 | 10 | 4 | 35 | 35 | 1 | 134 | 105 | 0.187 | 0.714 | 218 | 50290 | 28803 | 8.75 | 172 | 55647 | 38.4 | 36.7 | 0 / 0 / 0 | yes |
| 23 | 10 | 460 | 1 | 30 | 4 | 33 | 33 | 1 | 258 | 213 | 0.0814 | 0.636 | 218 | 1.45e+05 | 43388 | 8.25 | 142 | 1.55e+05 | 107 | 102 | 0 / 0 / 0 | yes |
| 23 | 10 | 460 | 2 | 3 | 4 | 29 | 29 | 1 | 100 | 78 | 0.2 | 0.69 | 233 | 38360 | 24093 | 7.25 | 105 | 41953 | 29 | 27.7 | 0 / 0 / 0 | yes |
| 23 | 10 | 460 | 2 | 10 | 4 | 38 | 38 | 1 | 191 | 158 | 0.152 | 0.763 | 216 | 2.83e+05 | 45453 | 9.5 | 171 | 2.88e+05 | 199 | 190 | 0 / 0 / 0 | yes |
| 23 | 10 | 460 | 2 | 30 | 4 | 33 | 33 | 1 | 270 | 225 | 0.0778 | 0.636 | 218 | 2.22e+06 | 45562 | 8.25 | 142 | 2.23e+06 | 1542 | 1473 | 0 / 0 / 0 | yes |
| 23 | 20 | 920 | 0 | — | 4 | 72 | 72 | 1 | 0 | 0 | — | 0 | 62.7 | 9660 | 0 | 18 | 501 | 12093 | 8.35 | 7.98 | 0 / 0 / 0 | yes |
| 23 | 20 | 920 | 1 | 3 | 4 | 65 | 65 | 1 | 124 | 88 | 0.282 | 0.538 | 141 | 66140 | 31905 | 16.2 | 428 | 72310 | 49.9 | 47.7 | 0 / 0 / 0 | yes |
| 23 | 20 | 920 | 1 | 10 | 4 | 61 | 61 | 1 | 198 | 165 | 0.131 | 0.426 | 128 | 1.97e+05 | 44795 | 15.2 | 305 | 2.06e+05 | 142 | 136 | 0 / 0 / 0 | yes |
| 23 | 20 | 920 | 1 | 30 | 4 | 64 | 64 | 1 | 306 | 256 | 0.049 | 0.234 | 102 | 5.72e+05 | 45561 | 16 | 385 | 5.89e+05 | 407 | 389 | 0 / 0 / 0 | yes |
| 23 | 20 | 920 | 2 | 3 | 4 | 69 | 69 | 1 | 155 | 113 | 0.265 | 0.594 | 147 | 1.5e+05 | 43266 | 17.2 | 442 | 1.57e+05 | 108 | 103 | 0 / 0 / 0 | yes |
| 23 | 20 | 920 | 2 | 10 | 4 | 61 | 61 | 1 | 200 | 167 | 0.13 | 0.426 | 128 | 1.12e+06 | 45562 | 15.2 | 305 | 1.13e+06 | 781 | 746 | 0 / 0 / 0 | yes |
| 23 | 20 | 920 | 2 | 30 | 4 | 64 | 64 | 1 | 306 | 256 | 0.049 | 0.234 | 102 | 8.87e+06 | 45562 | 16 | 385 | 8.88e+06 | 6136 | 5862 | 0 / 0 / 0 | yes |
| 31 | 5 | 310 | 0 | — | 8 | 45 | 44 | 1.02 | 0 | 0 | — | 0 | 213 | 930 | 0 | 5.5 | 77.5 | 2774 | 2.31 | 2.54 | 0 / 0 / 0 | yes |
| 31 | 5 | 310 | 1 | 3 | 8 | 44 | 44 | 1 | 103 | 73 | 0.262 | 0.614 | 267 | 5935 | 4239 | 5.5 | 73.1 | 9172 | 7.65 | 8.38 | 0 / 0 / 0 | yes |
| 31 | 5 | 310 | 1 | 10 | 8 | 43 | 43 | 1 | 162 | 121 | 0.185 | 0.698 | 277 | 17380 | 11320 | 5.38 | 67.9 | 21752 | 18.1 | 19.9 | 0 / 0 / 0 | yes |
| 31 | 5 | 310 | 1 | 30 | 8 | 43 | 43 | 1 | 293 | 241 | 0.102 | 0.698 | 275 | 50055 | 20163 | 5.38 | 63.5 | 57566 | 48 | 52.6 | 0 / 0 / 0 | yes |
| 31 | 5 | 310 | 2 | 3 | 8 | 41 | 41 | 1 | 115 | 83 | 0.278 | 0.78 | 296 | 13375 | 9093 | 5.12 | 74.6 | 16652 | 13.9 | 15.2 | 0 / 0 / 0 | yes |
| 31 | 5 | 310 | 2 | 10 | 8 | 40 | 40 | 1 | 229 | 194 | 0.127 | 0.725 | 289 | 96430 | 22802 | 5 | 64.1 | 1.01e+05 | 84 | 92.1 | 0 / 0 / 0 | yes |
| 31 | 5 | 310 | 2 | 30 | 8 | 43 | 43 | 1 | 329 | 276 | 0.0942 | 0.721 | 278 | 7.52e+05 | 23216 | 5.38 | 62.2 | 7.6e+05 | 633 | 694 | 0 / 0 / 0 | yes |
| 31 | 10 | 620 | 0 | — | 4 | 39 | 39 | 1 | 0 | 0 | — | 0 | 105 | 3410 | 0 | 9.75 | 188 | 5533 | 4.61 | 4.19 | 0 / 0 / 0 | yes |
| 31 | 10 | 620 | 1 | 3 | 4 | 39 | 39 | 1 | 74 | 57 | 0.216 | 0.41 | 158 | 22700 | 12730 | 9.75 | 186 | 27611 | 23 | 20.9 | 0 / 0 / 0 | yes |
| 31 | 10 | 620 | 1 | 10 | 4 | 35 | 35 | 1 | 125 | 104 | 0.12 | 0.429 | 171 | 67250 | 21657 | 8.75 | 182 | 74277 | 61.9 | 56.3 | 0 / 0 / 0 | yes |
| 31 | 10 | 620 | 1 | 30 | 4 | 28 | 28 | 1 | 112 | 103 | 0.0536 | 0.214 | 165 | 1.95e+05 | 23209 | 7 | 119 | 2.08e+05 | 173 | 157 | 0 / 0 / 0 | yes |
| 31 | 10 | 620 | 2 | 3 | 4 | 33 | 33 | 1 | 93 | 75 | 0.183 | 0.515 | 191 | 51530 | 20166 | 8.25 | 159 | 56474 | 47.1 | 42.8 | 0 / 0 / 0 | yes |
| 31 | 10 | 620 | 2 | 10 | 4 | 34 | 34 | 1 | 133 | 112 | 0.113 | 0.441 | 176 | 3.8e+05 | 23216 | 8.5 | 173 | 3.87e+05 | 323 | 294 | 0 / 0 / 0 | yes |
| 31 | 10 | 620 | 2 | 30 | 4 | 28 | 28 | 1 | 112 | 103 | 0.0536 | 0.214 | 165 | 2.99e+06 | 23216 | 7 | 119 | 3.01e+06 | 2506 | 2280 | 0 / 0 / 0 | yes |
| 31 | 20 | 1240 | 0 | — | 4 | 60 | 60 | 1 | 0 | 0 | — | 0 | 67.1 | 13020 | 0 | 15 | 291 | 15990 | 13.3 | 12.1 | 0 / 0 / 0 | yes |
| 31 | 20 | 1240 | 1 | 3 | 4 | 59 | 59 | 1 | 61 | 56 | 0.082 | 0.0847 | 79.7 | 88780 | 22207 | 14.8 | 290 | 96878 | 80.7 | 73.4 | 0 / 0 / 0 | yes |
| 31 | 20 | 1240 | 1 | 10 | 4 | 53 | 53 | 1 | 65 | 61 | 0.0308 | 0.0377 | 80.8 | 2.65e+05 | 23216 | 13.2 | 276 | 2.77e+05 | 231 | 210 | 0 / 0 / 0 | yes |
| 31 | 20 | 1240 | 1 | 30 | 4 | 60 | 60 | 1 | 80 | 78 | 0 | 0 | 68.4 | 7.68e+05 | 23216 | 15 | 291 | 7.93e+05 | 661 | 601 | 0 / 0 / 0 | yes |
| 31 | 20 | 1240 | 2 | 3 | 4 | 58 | 58 | 1 | 61 | 55 | 0.0984 | 0.103 | 83.5 | 2.02e+05 | 23206 | 14.5 | 289 | 2.1e+05 | 175 | 159 | 0 / 0 / 0 | yes |
| 31 | 20 | 1240 | 2 | 10 | 4 | 53 | 53 | 1 | 65 | 61 | 0.0308 | 0.0377 | 80.8 | 1.51e+06 | 23216 | 13.2 | 276 | 1.52e+06 | 1270 | 1155 | 0 / 0 / 0 | yes |
| 31 | 20 | 1240 | 2 | 30 | 4 | 60 | 60 | 1 | 80 | 78 | 0 | 0 | 68.4 | 1.19e+07 | 23216 | 15 | 291 | 1.2e+07 | 9978 | 9076 | 0 / 0 / 0 | yes |
| 41 | 10 | 820 | 0 | — | 4 | 237043 | 42 | 5644 | 0 | 0 | — | 0 | 4.63e+06 | 4510 | 0 | 10.5 | 200 | 5.25e+07 | 70.8 | 420 | 0 / 0 / 0 | yes |
| 41 | 10 | 820 | 1 | 3 | 4 | 60004 | 44 | 1364 | 176 | 71 | 0.193 | 0.773 | 1.12e+06 | 29550 | 24600 | 11 | 223 | 1.53e+07 | 20.6 | 122 | 0 / 0 / 0 | yes |
| 41 | 10 | 820 | 1 | 10 | 4 | 35382 | 43 | 823 | 333 | 148 | 0.111 | 0.86 | 6.75e+05 | 87720 | 82000 | 10.8 | 211 | 9.08e+06 | 12.2 | 72.7 | 0 / 0 / 0 | yes |
| 41 | 10 | 820 | 1 | 30 | 4 | 22569 | 41 | 550 | 631 | 294 | 0.0634 | 0.976 | 4.52e+05 | 2.54e+05 | 2.46e+05 | 10.2 | 242 | 6e+06 | 8.09 | 48 | 0 / 0 / 0 | yes |
| 41 | 10 | 820 | 2 | 3 | 4 | 26736 | 41 | 652 | 172 | 109 | 0.227 | 0.951 | 5.35e+05 | 67680 | 61410 | 10.2 | 416 | 6.86e+06 | 9.25 | 54.9 | 0 / 0 / 0 | yes |
| 41 | 10 | 820 | 2 | 10 | 4 | 6709 | 43 | 156 | 382 | 315 | 0.113 | 1 | 1.28e+05 | 5.02e+05 | 4.92e+05 | 10.8 | 435 | **2.21e+06** | **2.98** | **17.7** | 0 / 0 / 0 | yes |
| 41 | 10 | 820 | 2 | 30 | 4 | 1796 | 44 | 40.8 | 846 | 782 | 0.052 | 1 | 33832 | 3.96e+06 | 3.93e+06 | 11 | 582 | 4.43e+06 | 5.97 | 35.5 | 0 / 0 / 0 | yes |
| 41 | 40 (control) | 3280 | 0 | — | 4 | 16142 | 155 | 104 | 0 | 0 | — | 0 | 3.39e+05 | 67240 | 0 | 38.8 | 2038 | 1.43e+07 | 19.2 | 114 | 0 / 0 / 0 | yes |
| 41 | 110 (control) | 9020 | 0 | — | 4 | 1971 | 405 | 4.87 | 0 | 0 | — | 0 | 37300 | 5.01e+05 | 0 | 101 | 14428 | 4.59e+06 | 6.19 | 36.7 | 0 / 0 / 0 | yes |

Every one of the 51 cells recovered and verified its planted logarithm on
every target (260 targets); no partial was rejected by the exact eliminator,
no merge hit a limit, and no merged row failed the independent
counted-group replay.

### What the cells say

**n = 23 and n = 31 (r ≈ 2²¹ and 2^20.5) cannot discriminate.** At these
rungs `|F|³/6 ≥ r` already at 5 columns, so the unmodified three-summand
table decomposes nearly every trial (trials per complete relation 1.0–1.6).
The collector returns at most one complete relation per trial, as the
classic oracle does, so it cannot push that below 1; what the large primes
do here is *displace* direct hits with merged rows (merged share 0.2–0.9:
the loop meets a partial at an earlier `i` than the direct hit and returns
the merged row first) at a higher cost per relation (5 columns, n = 31: 267–
296 relation adds per relation against 213) and a much larger setup. The
`F ∪ D` table also saturates: its entries cap at `r / 2n` = 45,562 (n = 23)
and 23,216 (n = 31), the number of signed Frobenius orbits in the subgroup,
from `|D|/|F| = 10` upward — the "bounded" domain is then a large fraction of
the group and the table enumerates every orbit, which is the `|D| = r` limit
§7.4 reasoned about. On these two rungs every large-prime cell is worse end
to end than the plain oracle: 2.2–274× at 5 columns, 4.7–544× at 10, 6–750×
at 20 in one-target total, the factor growing as `|D|/|F|` (single) and as
`(|D|/|F|)²` (double) with the table.

**n = 41 (r = 2³⁹), 10 columns: the one cell with headroom.** The plain
oracle needs 5,644 trials per complete relation (4.63·10⁶ additions per
relation, 5.25·10⁷ GAE for one target, 420× the counted rho). Against that,
at equal base and columns:

| variant | trials / rel. | ÷ | rel. adds / rel. | ÷ | setup adds | one-target total GAE | ÷ | S / rho |
|---|--:|--:|--:|--:|--:|--:|--:|--:|
| LP = 0 (plain `m = 3`) | 5,644 | 1 | 4.63·10⁶ | 1 | 4,510 | 5.25·10⁷ | 1 | 420 |
| LP = 1, \|D\| = 3\|F\| | 1,364 | 4.1 | 1.12·10⁶ | 4.1 | 29,550 | 1.53·10⁷ | 3.4 | 122 |
| LP = 1, \|D\| = 10\|F\| | 823 | 6.9 | 6.75·10⁵ | 6.9 | 87,720 | 9.08·10⁶ | 5.8 | 72.7 |
| LP = 1, \|D\| = 30\|F\| | 550 | 10.3 | 4.52·10⁵ | 10.2 | 2.54·10⁵ | 6.00·10⁶ | 8.8 | 48 |
| LP = 2, \|D\| = 3\|F\| | 652 | 8.7 | 5.35·10⁵ | 8.7 | 67,680 | 6.86·10⁶ | 7.7 | 54.9 |
| LP = 2, \|D\| = 10\|F\| | 156 | 36 | 1.28·10⁵ | 36 | 5.02·10⁵ | **2.21·10⁶** | **23.8** | **17.7** |
| LP = 2, \|D\| = 30\|F\| | 40.8 | 138 | 3.38·10⁴ | 137 | 3.96·10⁶ | 4.43·10⁶ | 11.9 | 35.5 |
| control: plain `m = 3`, 40 columns (= \|F ∪ D\| at 3\|F\|) | 104 | — | 3.39·10⁵ | — | 67,240 | 1.43·10⁷ | 3.7 | 114 |
| control: plain `m = 3`, 110 columns (= \|F ∪ D\| at 10\|F\|) | 4.87 | — | 3.73·10⁴ | — | 5.01·10⁵ | 4.59·10⁶ | 11.4 | 36.7 |

Readings, all at one rung, 4 targets per cell, one-target cold totals:

* The single large prime cuts trials per relation by 4.1×, 6.9× and 10.3×
  at `|D|/|F|` = 3, 10 and 30 — inside §7.4's prediction of "÷ 3–10 with
  `|D| = 10|F|`". Relation-phase additions fall by the same factor (the
  probe loop is unchanged; the extra cost per trial is lookups and one
  canonicalisation, charged in GAE), and the setup grows with `|F|·|D|/2n`.
* The double large prime cuts trials per relation a further 2.1×, 5.3× and
  13.5× over the single at equal `|D|`: 8.7×, 36× and 138× against the
  plain oracle. Its graph closes at mean depth 2.6–4.0 with 5.6–12.6
  source partials per completed row; merge yield is 5–23 % of partials fed.
  Its table is `|F ∪ D|² / 4n` rather than `|F|·|D| / 2n`: at `|D| = 30|F|`
  the table (3.96·10⁶ adds) is 89 % of the one-target total, which is why
  the best cell is `|D| = 10|F|`, not 30.
* **Against the equal-envelope control** the single large prime is not a
  gain: at `|D| = 3|F|` and `10|F|` it costs 1.07× and 1.98× what the plain
  oracle costs on a base the size of `F ∪ D` (40 and 110 columns) with the
  same or a smaller table. The double large prime is: 2.1× below the
  40-column control at `|D| = 3|F|` and 2.08× below the 110-column control at
  `|D| = 10|F|`, with a table of the same size (4.92·10⁵ against 5.01·10⁵
  entries), 11 matrix rows instead of 101, and 820-point probe loops instead
  of 9,020-point ones. That is the only measured gain in this directory,
  and it is a constant on one rung.
* Nothing here approaches the reference: the best cell is 17.7× the counted
  rho at 39 bits (the plain 10-column oracle is 420×, the 110-column control
  36.7×). No exponent is fitted — one rung with headroom is not a fit.

### Not run

* n = 41 at 5 and 20 columns. The 10-column sweep took 50 minutes of wall
  time on a host at load average ≈ 70–155 for about 1.5 CPU-minutes of work
  (the `ic` process held ≈ 8 % of a core); the 5-column sweep (≈ 48,000
  trials per relation for the plain oracle, 8 targets) and the 20-column
  sweep (whose `LP = 2, |D| = 30|F|` cell builds a 3·10⁷-entry table four
  times) were started as part of `run.sh 41` and stopped after the 10-column
  report was written, before the 20-column `ic` process had produced
  anything; nothing of theirs is recorded. Resource limit on the shared host,
  not a finding; the sweep files are in `sweeps/` and `./run.sh 41` runs them.
* The research-visuals report (graph + PDF) the repository asks of every run
  is not in this directory; the table above and `results.json` are the record
  until it is produced.

## Verdict against §7.4

§7 item 4 of `LEAD_LADDER_METHODOLOGY_20261010.md` closed the *single*
large prime on the Frobenius-folded table with two statements: (i) with the
residual collector as the `|D| = r` limit, the collector beats the direct
table only below the break-even `|F| < 2.4·(K r)^{1/4}` (≈ 100 points at n =
23, where every rung's base is larger); (ii) a bounded `D` lowers the entry
rate and the collision rate by the same factor, so it cannot cross over where
`|D| = r` does not.

* **(i) stands, and the measurement lands on both sides of it.** At n = 23
  and 31 the break-even is 91 and 83 points and every base measured (230–
  1,240 points) is above it: every large-prime cell lost, by 2.2–740×. At n =
  41 the break-even is `2.4·r^{1/4}` ≈ 2,070 points (3,680 with `K = 10`)
  and the 820-point base is below it: the single large prime won 3.4–8.8×
  against the same-column plain oracle, inside the predicted ÷ 3–10. The
  formula's regime split is what the data shows.
* **(ii) does not carry to the double-large-prime graph, but the control
  takes most of what it gains.** The double graph's partials are found by the
  `F × D` and `D × D` tables, whose entry rate grows with `|D|` and `|D|²`,
  not by membership of a residual in `D`; its trials per relation fell 36×
  at `|D| = 10|F|` where the single fell 6.9×. End to end, however, the
  single large prime is at best level with (1.07×) and typically behind
  (1.98×) the plain oracle on a base the size of `F ∪ D`, and the double
  large prime is 2.1× ahead of it at equal table size. So the double graph
  changes the §7.4 conclusion from "no gain over the table" to "a 2× constant
  over spending the same table on more columns, at one rung, with the matrix
  10× smaller" — **engineering**, by the AGENTS.md §3 test: `S` fell, the
  ratio to the floor and to rho moved by a constant on one rung, no exponent
  changed. The break-even formula itself was not re-derived for the double
  graph; one rung with headroom (and two degenerate ones) does not support a
  fit, and the n = 41 cells at 5 and 20 columns that would show how the gain
  moves with `|F|` were not run.

What would move this: the two missing n = 41 cells and one more rung with
`|F|³/6 ≪ r` (n = 53, `r` ≈ 2⁵¹, or n = 47) under the same sweep, to see
whether the 2× over the control holds, grows or shrinks with `|F|` and `r`;
and a priced many-target variant (the collector's `D` columns are a
reusable relation-log asset only if the base's logs are kept), reported
separately from the one-target figures above as the rules require.

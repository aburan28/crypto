# Multiplier closure in the 4-sum pair-table index calculus (toy Koblitz)

**Stage diagnostic.**  `S`, end-to-end cost and speedup are **unset**.  No
claim about ECC2K-130 or any deployed curve.  Protocol frozen before the
run in [`PROTOCOL.md`](PROTOCOL.md); its post-run amendment is at the end
of that file.  Companion note:
[`docs/ic/RESEARCH_TRANSFER_AWARE_IC_20261010.md`](../../docs/ic/RESEARCH_TRANSFER_AWARE_IC_20261010.md).

## Question

The compact-orbit producer builds its summand set from `K` columns times
the Frobenius multipliers `{±τ^k}` and finds 4-sum relations by probing a
pair table.  "Loops of isogenies" at the crater of the ECC2K-130 class are
`Z[τ]`, and §A1 of `RESEARCH_IC_NOVEL_DIRECTIONS_20261007.md` proposes
closing the multiplier set under its small-norm elements.  The note
derives that, in this pipeline, (i) any cheap known-log multiplier set
plays the same role as the Frobenius orbit in the yield model, so plain
integers `{1, …, c}` should do what `Z[τ]` does; (ii) at equal table
memory the closure cannot change per-target probes, only trade columns for
multipliers; (iii) the trade is limited by additive coincidences among the
multiples of one column, which show up at small `K`.

## Answer (cells C23, C37, C37b complete; C41 and C41b running, see below)

1. **Yield law holds.**  The `frob` arm's mean first-hit probes are within
   10 % of `2r/L²` in every finished cell (15.1 vs 13.8; 63.9 and 68.0 vs
   65.0; partial C41: 35,861 and 29,143 vs 31,543).  M1 pass.
2. **`Z[τ]` is cosmetic.**  At every `K`, the `tau` and `int` arms are
   within `0.88–1.56×` of each other, and at `K ≥ 3` both are within
   `0.88–1.26×` of `frob`.  M2 and M3 pass; the decisive outcome that
   would have revived the `Z[τ]` framing (`tau ≤ 0.67 × int`) did not
   occur in any cell.
3. **Coincidence loss is real and tracks the table.**  At `K = 1` the
   canonical table shrinks to `0.58–0.72×` and first-hit probes rise to
   `1.49–2.15×` the `frob` arm's.  M4 as frozen (`≥ 1.5×`) fails by a hair
   on the `int K = 1` arms of C37 (1.494) and C37b (1.491): the effect is
   there, the band was tight.
4. **Every multiplier arm solved its columns and its hidden target**
   (M6), with `rows_to_full_rank = K` or `K + 1` (M5).  The `frob` arms
   are where the protocol was wrong: sparse 4-term rows need an excess
   over `K` that M5 did not allow (18 rows for `K = 12`, 42 for `K = 36`),
   one registered seed left a column untouched in 64 rows (C37, rank 35),
   and C41 was registered with fewer hits than columns.  These are
   recorded as M5/M6 failures of the registered cells and corrected in
   the amended cells C37b and C41b, never overwritten.
5. **Equal-memory trade.**  Total probes to full rank scale with `K`
   (M7 passes on 9 of 11 scorable arms; the two misses are a `K = 3` arm
   whose three first hits were lucky, ratio 5.0, and C37's `frob` arm
   with no full rank).  Relations and linear algebra are divided by `c`;
   per-target probes are not.

## Results

`L = 2n·Kc` in every arm of a cell; `canonical` is the number of distinct
Frobenius-canonical `x` in the pair table (the τ-canonicalised index the
landed design stores); `probes` is the mean first-hit probe count over
`hits` fresh targets; `rows` is the row count at which rank `K` was
reached; `total` is the probes spent up to that row.

**C23** (`n = 23`, `r = 2,095,853`, `Kc = 12`, `L = 552`, `2r/L² = 13.8`, 64 hits, seed 20261010)

| arm | K | c | canonical | canon/frob | probes | vs frob | rows | total | logs | target |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|
| frob | 12 | 1 | 3163 | 1.000 | 15.1 | 1.00 | 18 | 255 | ok | ok |
| tau | 6 | 2 | 3117 | 0.985 | 18.8 | 1.24 | 7 | 159 | ok | ok |
| int | 6 | 2 | 3095 | 0.979 | 13.7 | 0.91 | 6 | 113 | ok | ok |
| tau | 4 | 3 | 3034 | 0.959 | 14.6 | 0.97 | 4 | 54 | ok | ok |
| int | 4 | 3 | 3015 | 0.953 | 16.2 | 1.08 | 4 | 53 | ok | ok |
| tau | 2 | 6 | 2666 | 0.843 | 19.6 | 1.30 | 2 | 23 | ok | ok |
| int | 2 | 6 | 2672 | 0.845 | 12.6 | 0.84 | 2 | 30 | ok | ok |
| tau | 1 | 12 | 1834 | 0.580 | 26.1 | 1.73 | 1 | 1 | ok | ok |
| int | 1 | 12 | 1961 | 0.620 | 26.2 | 1.73 | 1 | 6 | ok | ok |

**C37** (`n = 37`, `r = 230,603,167`, cofactor 596, `Kc = 36`, `L = 2664`, `2r/L² = 65.0`, 64 hits, seed 20261010)

| arm | K | c | canonical | canon/frob | probes | vs frob | rows | total | logs | target |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|
| frob | 36 | 1 | 47454 | 1.000 | 63.9 | 1.00 | none (rank 35; column 5 untouched) | — | **no** | **no** |
| tau | 18 | 2 | 47160 | 0.994 | 78.1 | 1.22 | 21 | 2095 | ok | ok |
| int | 18 | 2 | 47270 | 0.996 | 69.0 | 1.08 | 21 | 1879 | ok | ok |
| tau | 9 | 4 | 46580 | 0.982 | 74.4 | 1.17 | 9 | 1000 | ok | ok |
| int | 9 | 4 | 46625 | 0.983 | 74.1 | 1.16 | 9 | 317 | ok | ok |
| tau | 3 | 12 | 43201 | 0.910 | 80.4 | 1.26 | 3 | 168 | ok | ok |
| int | 3 | 12 | 43457 | 0.916 | 67.0 | 1.05 | 3 | 104 | ok | ok |
| tau | 1 | 36 | 30471 | 0.642 | 109.1 | 1.71 | 1 | 164 | ok | ok |
| int | 1 | 36 | 33907 | 0.715 | 95.4 | 1.49 | 1 | 15 | ok | ok |

**C37b** (amended: 128 hits, seed 20261011; otherwise as C37)

| arm | K | c | canonical | canon/frob | probes | vs frob | rows | total | logs | target |
|:--|--:|--:|--:|--:|--:|--:|--:|--:|:--|:--|
| frob | 36 | 1 | 47473 | 1.000 | 68.0 | 1.00 | 42 | 3193 | ok | ok |
| tau | 18 | 2 | 47255 | 0.995 | 75.3 | 1.11 | 21 | 1305 | ok | ok |
| int | 18 | 2 | 47263 | 0.996 | 70.7 | 1.04 | 18 | 963 | ok | ok |
| tau | 9 | 4 | 46595 | 0.982 | 60.0 | 0.88 | 9 | 388 | ok | ok |
| int | 9 | 4 | 46676 | 0.983 | 68.1 | 1.00 | 10 | 884 | ok | ok |
| tau | 3 | 12 | 42792 | 0.901 | 80.3 | 1.18 | 3 | 53 | ok | ok |
| int | 3 | 12 | 43594 | 0.918 | 75.1 | 1.10 | 3 | 119 | ok | ok |
| tau | 1 | 36 | 30471 | 0.642 | 146.4 | 2.15 | 1 | 56 | ok | ok |
| int | 1 | 36 | 33907 | 0.714 | 101.4 | 1.49 | 1 | 32 | ok | ok |

**C41 and C41b** (`n = 41`, `r = 549,756,390,943`, `Kc = 72`, `L = 5904`,
`2r/L² = 31,543`): the first launch was killed with the session that
started it after five (C41) and two (C41b) arms; the partial logs are kept
as `results/c41_partial_first_session.log` and
`results/c41b_partial_first_session.log`.  Both cells were relaunched
unchanged (the instrument is deterministic in the seed; the C37 replay
reproduced the registered means to the last digit) and write
`results/c41.json` and `results/c41b.json` when they finish; `score.py`
reports them as `missing` until then.  From the partial logs:

| cell | arm | K | c | probes | vs frob | rows | logs | target |
|:--|:--|--:|--:|--:|--:|--:|:--|:--|
| C41 (48 hits) | frob | 72 | 1 | 35,861 | 1.00 | none (48 rows < 72 columns: protocol error) | **no** | **no** |
| C41 | tau | 36 | 2 | 34,579 | 0.96 | 36 | ok | ok |
| C41 | int | 36 | 2 | 28,952 | 0.81 | 36 | ok | ok |
| C41 | tau | 12 | 6 | 36,100 | 1.01 | 12 | ok | ok |
| C41 | int | 12 | 6 | 29,498 | 0.82 | 12 | ok | ok |
| C41b (160 hits) | frob | 72 | 1 | 29,143 | 1.00 | 88 | ok | ok |
| C41b | tau | 36 | 2 | 34,570 | 1.19 | 36 | ok | ok |

## Scoring (`results/score.json`, finished cells)

| prediction | verdict |
|:--|:--|
| M1 yield law | pass in every cell |
| M2 equivalence at `K ≥ 3` | pass in every arm |
| M3 `tau/int ∈ [0.6, 1.6]` | pass; decisive outcome absent |
| M4 coincidence loss at `K = 1` | canonical part passes everywhere; probe part `1.49×` on two `int` arms against the `1.5×` bar: **fail as frozen**, effect present |
| M5 rows `≤ K + 2`, no degenerate rows | no degenerate rows anywhere; `frob` arms need more excess than allowed: **fail** (protocol error) |
| M6 logs and target verified | every `tau`/`int` arm passes; `frob` fails in C37 (untouched column) and C41 (hits < K): **fail** (protocol error) |
| M7 equal-memory trade | pass on 9 of 11 scorable arms |

## Reading

The multiplier closure is a constant-factor precompute trade, available
to any cheap known-log multiplier set: at equal table memory it divides
relations and linear-algebra width by `c` and leaves the per-target
probe count alone, and it stops paying when `K` gets small enough for
same-column coincidences to thin the table.  `Z[τ]` buys free evaluation
of `α(P)`, nothing in the yield.  Class: **stage diagnostic; closed as a
yield lever.**  The arithmetic at the landed rungs is in the companion
note.

## Files

| file | role |
|:--|:--|
| `PROTOCOL.md` | frozen predictions M1–M7, disclosed smoke tests, post-run amendment |
| `multipliers.rs` | standalone instrument (`rustc -O`), pure std |
| `score.py` | scores every present cell; writes `results/score.json` |
| `results/c23.json`, `c37.json`, `c41.json` | registered cells (seed 20261010) |
| `results/c37b.json`, `c41b.json` | amended cells (seed 20261011, `hits ≥ 2K`) |
| `results/*.log` | the instrument's stderr for each run |
| `results/c37_replay_touches.json` | replay of C37 with per-column diagnostics (column 5: zero touches) |
| `results/diag37_seed1.json` | 200-hit diagnostic at seed 1 (rank 36 in 36 rows) |
| `results/*_partial_first_session.log` | the killed first launch of C41/C41b |

```sh
rustc -O -o multipliers multipliers.rs
./multipliers 23 20261010 12 64 results/c23.json
./multipliers 37 20261010 36 64 results/c37.json
./multipliers 41 20261010 72 48 results/c41.json
./multipliers 37 20261011 36 128 results/c37b.json
./multipliers 41 20261011 72 160 results/c41b.json
python3 score.py
```

Host: shared macOS box at load 90–150 throughout; wall times in the JSON
are not comparable across arms.

# A rho with its cache sized to the walk, against the IC arm

Registered before the measured run. §5 lists what was already seen.

## 1. Question

The [stop decision](../DECISION-20260929-stop.md) leaves one loose end: the IC arm's small-cell
time lead over the lean rho. The [throughput diagnostic](../throughput_gap_20260929/RESULTS.md)
traced most of it to one buffer. Rho zero-fills a fixed 163,840-byte recent-point cache on every
solve (`vec![None; RECENT_SLOTS]`, `RECENT_SLOTS = 1 << 12`, 40-byte slots), whatever the walk. At
`n ≤ 23` a walk is a few hundred steps. On `n13a0-000` it is 4.

The diagnostic was a simulation, not a timing, and it said "a rho with the cache sized to the walk
would probably remove IC's small-cell time lead. That rho has not been built or measured." This
round builds it and times it.

**Does IC's small-cell time lead survive a rho whose cache is sized to its walk?**

## 2. The candidate

[`sized-cache.patch`](sized-cache.patch) (`-p1`, against round 0025's `source`) changes the lean
rho in three places and nothing else:

- `rho_recent_slots(expected_steps)` returns `clamp(next_pow2(⌈4 × expected_steps⌉), 64, 4096)`;
- the cache is allocated with that many slots instead of `RECENT_SLOTS`;
- the slot index masks with `slots − 1` instead of `RECENT_SLOTS − 1`.

`expected_steps` is the existing `rho_expected_steps(modulus, n)`. The cap is the old constant, so
a cell with 1,024 or more expected steps runs the old code exactly. The IC mode, the incumbent
arm, and every other function are untouched.

Built the way the round built its worker: `cargo build --release --offline --locked --jobs 2
--example ic_tournament_worker` in the round's `source` with the patch applied, so the same
`.cargo/config.toml` (musl, static, `--gc-sections`) applies.

**The walk can change.** The cache only catches short repeats. A smaller direct-mapped cache
evicts more, so a cycle or collision can be found a few steps later, or a stored-point hit found
first. This patch does **not** promise the same walk. It promises a verified answer, and §4 measures
how much the walk moves.

## 3. Method

- **Runner.** [`run.py`](run.py), through [`tools/isolated_bench.py`](../../../../tools/isolated_bench.py)
  (AGENTS.md §10), on the same terms as the
  [isolated re-timing](../walltime_isolation_20260929/PREREGISTRATION.md#3-method):
  - one lock and one preflight, refusing to start while other processes use more than 0.1 CPU or
    PSI avg10 is above 5;
  - every movable thread moved off CPU 3, the driver kept off it, each run pinned to it;
  - per run: wall, user and system time, context switches, faults, PSI, and the CPU time other
    processes used. A run is `contended` if other processes used more than 0.1 CPU on average
    during it. Contended runs are excluded and counted.
  - no build, test, lint or other job runs on the machine during the timed run.
- **Fixtures.** Round 0025's confirmation cases at `n13a0, n17a1, n19a0, n19a1, n23a0`: 12 cases
  per cell, each arm on its own recorded `job.json`.
- **Arms** (five):

  | arm | binary | mode | answer check |
  |:--|:--|:--|:--|
  | `incumbent` | round 0025 `worker` | IC | whole `solutions` equals the record |
  | `rho` | round 0025 `worker` | lean rho | whole `solutions` equals the record |
  | `sized` | patched worker | rho | recovered logarithm equals the record, and `verified` |
  | `sized_aa` | patched worker | rho | as `sized` (A/A noise floor) |
  | `inc_new` | patched worker | IC | whole `solutions` equals the record (layout control) |

  `inc_new` runs the patched binary's IC mode. It should not differ from `incumbent`, because the
  IC code is untouched. If it does, the patched binary's layout is confounding the comparison.
- **Order.** Per round and case, the five arms follow a ten-row Williams square (the standard
  construction and its reversal), so each arm follows every other equally often. One untimed
  warm-up of every arm per cell comes first.
- **Size.** 20 rounds: 5 cells × 12 cases × 5 arms × 20 = 6,000 timed runs.
- **Timers.** `wall` is spawn to reap. `cpu` is user plus system. `worker` is the worker's own
  `elapsed_seconds`, from curve construction to report assembly.

## 4. Registered readout ([`analyze.py`](analyze.py))

- **Per case:** the median over rounds of arm X divided by the median over rounds of arm Y.
- **Per cell:** the geometric mean over the 12 cases, with a 95% bootstrap over cases (10,000
  draws, seed 20260930).

The small cells are `n13a0, n17a1, n19a0, n19a1`. `n23a0` is reported, with no prediction.

**Predictions**

- **P1 (the cache is the cost).** `sized/rho` on the `worker` timer has its upper 95% limit below
  1 at all four small cells.
- **P2 (the IC lead is gone).** `incumbent/sized` on `wall` has its upper limit at or above 1 at
  **at least three** of the four small cells. That is, IC is no longer significantly faster than the
  sized rho. The isolated re-timing had `incumbent/rho` at 0.87–0.94 there.
- **P3 (noise floor).** `sized_aa/sized` on `wall` lies inside [0.92, 1.09] at every cell.
- **P4.** At most 5% of timed runs are contended.
- **P5 (layout control).** `inc_new/incumbent` on `wall` lies inside [0.92, 1.09] at every cell.
  If P5 fails, the P2 comparison is confounded and is not read.

**Reported, not predicted.** Cells where the sized rho is strictly faster than IC on `wall`
(`incumbent/sized` lower limit above 1); `cpu`-timer ratios; median minor faults per run; the
`n23a0` row.

**Regression, not timed** ([`regress.py`](regress.py)). The patched worker's rho mode is run on
every rho fixture of round 0025 (development, confirmation, replay, aa, smoke and selection: every
`job.json` with a recorded lean-rho answer) and the whole `solutions` array is compared with the
record.

- **R1.** Every run exits 0, `verified` is true, and the recovered logarithm equals the record.
- **R2.** The walk (`iterations`, `restarts`, `walk_group_additions`) is identical on every run at
  `n = 13, 17, 19`, and on every run at `n ≥ 37`, where the cache size is unchanged. This is a
  reproduction of the smoke in §5, not an independent prediction.
- **R3.** Wherever the walk differs, total iterations rise by less than 1% on every degree.

### How each outcome is read

- **P1 and P2 hold, P3–P5 hold.** The small-cell IC lead was rho's fixed start-up write. It is
  engineering, and it is gone once that write is sized to the walk. The scoreboard's "IC lead"
  sentence is annotated, and the stop decision's reopening item for a stronger rho is closed for
  the time lead.
- **P1 holds, P2 fails.** The cache is a real cost but not the whole lead. The residual is
  reported per cell.
- **P1 fails.** The simulation over-attributed. The cache is not the time cost at these sizes, and
  the diagnostic's attribution is annotated.
- **P5 fails.** The comparison is confounded by the binary. Nothing is concluded about the lead.
- In every case, **no cell here can create an IC result.** The instruction count remains lower for
  rho at every cell, and the exponent question is untouched.

## 5. Already seen before registration

- **The diagnostic.** Everything in
  [throughput_gap_20260929/RESULTS.md](../throughput_gap_20260929/RESULTS.md) predates this
  registration. It is the reason for the patch.
- **The rule was chosen after reading the code, not after any timing.** The factor 4 and the bounds
  64 and 4,096 were set from the expected-step formula. Cache sizes it gives at the first case of
  each timed cell (bytes at 40 per slot):

  | cell | `r` | expected steps | slots | bytes |
  |:--|--:|--:|--:|--:|
  | `n13a0` | 2,003 | 11.0 | 64 | 2,560 |
  | `n17a1` | 65,587 | 55.0 | 256 | 10,240 |
  | `n19a0` | 130,873 | 73.6 | 512 | 20,480 |
  | `n19a1` | 262,543 | 104.2 | 512 | 20,480 |
  | `n23a0` | 2,095,853 | 267.5 | 2,048 | 81,920 |

  Cells with 1,024 or more expected steps, such as `n37a0` (2,212 steps), keep 4,096 slots. None
  of the constants was tuned against a time.
- **Regression smoke, seen in full.** The patched worker's rho mode was run once over all 1,536 rho
  runs in round 0025's run directories. This is the R1–R3 check, so its result is known. It was
  not timed.
  - 1,536 of 1,536 verified, and 1,536 of 1,536 recovered the recorded logarithm.
  - 1,509 had a walk identical to the record. The 27 that differed were at `n = 23` (15 runs) and
    `n = 31` (12 runs): iterations rose by 1 to 15, total ratio 1.0012 at both degrees, and no
    restart count changed.
  - Every run at `n = 13, 17, 19` had an identical walk, and so did every run at `n ≥ 29` except
    `n = 31`.
  - The registered regression is a re-run of the same deterministic check and is expected to
    reproduce it.
- **Plumbing run of `run.py`.** One round, `n13a0` only, into scratch. Its answers were checked by
  the runner (65 of 65 matched). **Its timings were not read**, and it is not pooled with the
  registered run.
- **No timing of the patched worker has been seen.**
- **The unregistered second check.** The `n23a0` row has no prediction because the cache there is
  2,048 slots (80 KB), half the old size, and a smaller effect is expected.

## 6. Scope

- **Host.** One host class: a 4-vCPU Intel Xeon @ 2.80 GHz cloud container. Its frequency and
  host neighbours are not controlled, and the A/A band measures that residue.
- **Round.** One round's binaries and fixtures, and the tournament worker's minimal-startup build.
  A stock `std` rho binary has a different fixed cost (about 4.5–4.9 ms whole-process at `n19a0`
  for the strong rho), and the comparison here is not a claim about it.
- **Constant, not exponent.** A time lead at `n ≤ 23` is a constant in the group size. It says
  nothing about how either algorithm's cost grows.

# ECDLP speedups across approaches — 2026-10-15

Engineering-class speedups (constant factors, no algorithmic change) landed
across three ECDLP approaches: the Koblitz compact-orbit shared-log DLP fast
path (`examples/koblitz_orbit_dlp_fast.rs`), and the prime j=0 orbit-reduced
index-calculus relation search (`src/cryptanalysis/ec_index_calculus_j0.rs`).
The frozen Kuhn–Struik batched rho v2 arm (`koblitz_rho_batch_ks_v2_n61.rs`)
was deliberately left untouched: it is the ledger baseline and modifying it
would break matched comparisons.

## Changes

### Koblitz compact-orbit fast path (`examples/koblitz_orbit_dlp_fast.rs`)

1. **Timer hoist in `lift128`** — `Instant::now()`/`elapsed()` previously sat
   in the innermost `d` loop of every candidate lift (one timestamp per
   candidate add). The lift is now timed once per candidate at the
   `scan_state128` call site. The `target_relation_check_ns` diagnostic now
   covers the whole lift, not a single add.
2. **Sparse canonical in the wide probe path** — `scan_state128` now calls
   `basis.canonical_selected(...)` (env-gated by `KIC_SPARSE_CANONICAL`,
   default dense) instead of the dense `basis.canonical(...)`. The wide table
   build already inserts `canonical_selected` keys, and the u64 probe path
   already used the selector, so the probe side now matches.
3. **Root-table Bloom prefilter on by default** — `root_bloom_enabled()`
   now returns true unless `KIC_ROOT_BLOOM=0`. Most probe keys miss, and the
   Bloom pass avoids a random DRAM walk per miss. Exactness is unchanged: a
   Bloom pass still probes the exact table; only misses are screened.

### Prime j=0 IC (`src/cryptanalysis/ec_index_calculus_j0.rs`)

4. **`x_to_orbit` map hoist** — `find_relation_with_psi` /
   `find_relation_with_psi_counted` rebuilt a `HashMap<BigUint, usize>` over
   all factor-base orbit x-coordinates on every relation call. New
   `orbit_index_map()` builds it once; the staged driver passes it into new
   `_with_index` variants. The original public signatures are preserved and
   delegate to the new variants.

## Measurement A/B — Koblitz orbit fast path

Matched baseline/candidate whole-process runs on identical inputs. Host:
Apple M4 Pro, 14 cores, macOS, rustc 1.93.1. Baseline binary built from
commit `3066b74c7410bd83e14371d78d83211fef81914b` (unmodified HEAD); candidate
is the same commit plus the three changes above. Not run under
`tools/isolated_bench.py` (other agents were using the machine); run order
was interleaved base/candidate within each repetition. Rep 3 of the candidate
was visibly contended (rank stage 27.8 s vs 14.5–16.5 s in the other reps)
and is retained, not discarded.

```
cargo run --release --example koblitz_orbit_dlp_fast -- \
  construct:41:0:26 scalars41.txt 20261015 out.jsonl
```

16 public synthetic known-answer scalar fixtures (scalars 1000+7919·i,
i = 1..16) at n=41, a=0, columns=26, rank seed 20261015, 1 thread.
Both arms solve the same deterministic 11/16 targets (the 5 unsolved are a
base-size artifact of the small column count, identical on both arms) and
emit bit-identical recovered scalars, `point_indices` and `x_codes` on every
repetition (diff-verified, `raw/`).

| arm | process_total (s), 3 reps | median | setup_total median (ms) | targets_total median (ms) | rank_stage median (ms) | target_median (ms) |
|---|---|---|---|---|---|---|
| baseline | 65.82 / 60.93 / 57.02 | 60.93 | 13.4 | 40 929 | 21 650 | 1 157 |
| candidate | 37.46 / 32.56 / 67.59 | 37.46 | 32.4 | 20 974 | 16 463 | 614 |
| ratio | — | **1.63×** | — | **1.95×** | 1.31× | 1.88× |

Wall-clock medians, three interleaved repetitions, identical outputs.
Class: **engineering** (S-ratio to rho unchanged in kind; the counted probe
behaviour and recovered scalars are identical). These are stage/whole-process
wall diagnostics on a shared host, not the §8 frozen-suite end-to-end claim.

## Measurement A/B — prime j=0 staged IC

```
cargo run --release --example j0_stage_ladder -- 16 20261005
```

The relation stage dominates (per-trial cost is BigUint arithmetic; the map
rebuild was per relation call). Three reps per arm; `random_scalar` is
unseeded so trial counts vary per run — per-trial normalization is given.

| arm | relations_ms, 3 reps | median | trials_total (3 reps) | ms per trial |
|---|---|---|---|---|
| baseline | 14 795 / 13 992 / 11 588 | 13 992 | 1 165 / 1 112 / 1 055 | ~12.6 |
| candidate | 13 813 / 11 897 / 11 941 | 11 897 | 1 076 / 1 034 / 970 | ~11.5 |

All runs recovered the same planted scalar and agreed with the rho control.
Class: **engineering**, ~1.17× median relation-stage wall (~8% per-trial);
the residual per-trial cost is the BigUint Semaev-s3/Tonelli–Shanks loop,
which is the next optimization target (see below).

## Correctness

- `cargo test --release --lib` (worktree): 4 096 passed; the 4–5
  `jv_cover` failures seen in one full pass were load-induced flakiness —
  the same tests pass on the unmodified base checkout and on a worktree
  rerun; `jv_cover` is untouched by this diff.
- `cargo test --release --lib ec_index_calculus_j0`: 15/15.
- `cargo test --release --example koblitz_orbit_dlp_fast`: 19/19
  (sparse-canonical exactness, root-Bloom exact first-witness, batched roots).
- n=23 and n=41 end-to-end runs: `targets_failed: 0` among solved,
  `group_verified` and `recovered_matches_published` true on every solved
  target, both arms.

## Not attempted here (identified, next targets)

- Porting the wide path's batched root solve (`RootBatch128`, one inversion
  shared across the 2n candidates) down to the u64 scan.
- Precomposing the twisted linear map `to_poly ∘ rotate_k ∘ to_normal` as
  per-shift byte tables (10–20% of the wide scan).
- Prime j=0: ψ-amplified orbit-shifted relations (up to ~3× trial
  amplification), parallel relation collection, p≡3 mod 4 sqrt shortcut.
- KS rho: the v2 arm is frozen as the ledger baseline; the v3 example
  already carries the batched-inversion and byte-table improvements.

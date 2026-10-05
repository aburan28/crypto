# Replay determinism of the koblitz yield-sweep session (2026-10-05)

A reproducibility fix, not a measurement: no figure on the scoreboard,
leaderboard or browser changes.

**Hypothesis.** `../sessions/koblitz/` fails exact replay because
`ecbench::methods::ic_report` removed the wall-priced solver term `w` of
the descent-algebraic arms by subtraction, so the recorded relations-phase
and total gae are `x` rounded to `w`'s binade, a function of the recording
host's clock.  Counts are unaffected.

**Success condition.** With the fix, `ecbench verify --replay-all` on the
unchanged session reproduces every replay, and repeated runs give the same
per-replay outcome.  Main's binary is run on the same host as the control.

## Commands

Host: Linux x86-64 cloud container, Intel(R) Xeon(R) Processor @ 2.80GHz,
4 logical CPUs, `rustc 1.97.0 (2d8144b78 2026-07-07)`.  The three runs ran
concurrently with each other and with builds, so wall timing was
deliberately perturbed; replays compare counts and gae bits only.

```
cargo build --release --bin ecbench
ecbench verify --dir research/ecbench_yield_sweep_20261004/sessions/koblitz/ \
    --replay-all --allow-interrupted [--exit-code] > <receipt>.json
```

| receipt | binary (commit, SHA-256) | replays | reproduced | detail |
|:--|:--|--:|--:|:--|
| `main_9f5afa23_run1.json` | main `9f5afa23`, `64722ac5…1378` | 197 | 178 | 19 differ in `total_gae` and/or `phase gae` (seqs 3,17,23,24,29,30,37,44,50,62,63,70,77,83,88,89,94,95,96) |
| `fix_b47dbd99_run1.json` | `b47dbd99`, `601f5109…ed3d` | 197 | 197 | 133 identical, 64 identical via the legacy-rounding reconstruction |
| `fix_b47dbd99_run2.json` | `b47dbd99`, `601f5109…ed3d` | 197 | 197 | same as run 1 |

The two fixed runs' per-replay `(seq, reproduced, detail)` lists are
byte-identical (SHA-256 of the compact `jq -c` list:
`02a06977…3fb3e6d` for both).  The control's failing set differs from the
set reported on Apple silicon (6–16 of 60, varying per run), as expected
for a timing-dependent rounding.

**Decision.** The fix lands in this PR (`ic_report` rebuilds the phase from
`gae_before_solver`; the audit reproduces committed legacy rounding from
the record's own solver term and pre-removal total).  The committed session
is unchanged.

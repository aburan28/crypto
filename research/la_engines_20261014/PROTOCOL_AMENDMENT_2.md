# Amendment 2 (2026-10-14, before any data for these arms)

Recorded while amendment 1 was running; it adds arms and changes none of
the frozen ones.

## Why

1. **CADO balancing abort.** At D = 32,768 seed 3, CADO `bwc`'s internal
   balancing aborted on the `replay-dl` matrix ("The current code expects a
   column permutation replicated on rows, not the converse … pass
   `--reorder columns` to mf_bal"). It is matrix-dependent (the `merge-dl -t
   4` output varies run to run) and reproduced with one thread. Under the
   amendment-1 stop rule the failure removed `cado-full` at every D ≥ 32,768.
   Rebalancing with `mf_bal --reorder columns` (what the message asks for)
   and rerunning `bwc` at the same 2×2 split solved and verified it
   (diagnostic runs, not data; logged in `raw/contention.log`).
2. **Our merge thresholds.** Amendment-1 data at D ≤ 16,384 show CADO's
   `merge-dl` (default `target_density` 100) leaves a far smaller matrix
   than our filter (defaults: merge columns of weight ≤ 8, rows ≤ 32
   nonzeros). Whether more aggressive merging helps our pipeline is the
   cheapest lever this study can test.

## Added arms (m = 3, `--coeffs positive`, seeds 1–3, four threads,
1,800 s timeout, same stop rule)

| arm | what | D |
|:--|:--|:--|
| `cado-full-fb` | amendment-1 `cado-full`, plus: on that specific balancing abort, `mf_bal --reorder columns` into a fresh wdir and rerun `bwc` at 2×2; `reorder_retry` recorded; the timed `bwc` step is `mf_bal` + rerun | 1,024 … 262,144 |
| `full-bw-m32r100` | `full-bw` with `merge_max_weight` 32, `max_row_weight` 100 | 16,384 … 262,144 |
| `full-bw-m128r200` | `full-bw` with `merge_max_weight` 128, `max_row_weight` 200 | 16,384 … 262,144 |

The pipeline comparison uses `cado-full-fb`; amendment-1 `cado-full` rows
stay in the record and the report says where they agree.

Binaries: `cado-full-fb` uses the `/tmp/cbt` build of `abdf48826` (the
amendment-1 instrument); the tuning arms use the `/tmp/cbt2` build of the
commit that adds `--merge-max-weight` / `--max-row-weight` (no other change).

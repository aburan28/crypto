# F1 on p₂, second run: annex

F1 (Amendment 2) was run twice on p₂ = 2013265921 with the same seed rule:

- **A**, the F1 rows of record: `runs/F1-{A0,A11}-m{2,3}-p2.jsonl` (parent
  lane, commits `f287673fe` and `fa60aba73`).
- **B**, run-20261015T2307Z (experiment lane, harness at `d5f675260`,
  2026-10-15 23:07Z to 2026-10-16 about 00:08Z, Xeon @ 2.10 GHz, 4 cores).
  Its report is `reports/run-20261015T2307Z/`.

Run B's full rows were committed at `runs/p2/` in `cb750c828`
(`git show cb750c828:research/factor_base_designs_20261015/runs/p2/A11-m2.jsonl`).
The F1 rows of record already contain every system B measured, so B's
duplicate rows were removed when the lanes were reconciled. What this
directory keeps:

| file | what it is |
|:--|:--|
| `crosscheck.txt` | `fb_designs_slope compare-rows A B` on all four series, run before B's rows were removed |
| `A11-m2-p2-t10-t11.jsonl` | B's 4 rows that A does not hold as completed: A11, `m = 2`, `t = 10` (both `ok`, D = 4) and `t = 11` (both timeouts at 300 s, partial D = 4) |
| `*-p2.log` | B's run logs, one line per system |

**What the cross-check says.**
- 52 systems are shared. 50 are identical in construction (seeds, `x_R`,
  curve, planted abscissae, counts) and F4 output (D, `degree_reached`,
  `max_cols`, `inconsistent`, `planted_ok`).
- The other 2 are A11, `m = 2`, `t = 10`. Their construction is identical,
  but A timed out at 300 s (305 s and 313 s, partial D = 4) and B completed
  (259 s and 256 s, D = 4). That is a host-speed difference, not a
  disagreement.
- B continued that series to `t = 11`, which timed out.

**What it does not change.** F1's decision rule is registered for `m = 3`, and
every `m = 3` row is identical in A and B. The `m = 2` series is descriptive.
B's extra completed cell (`t = 10`) would move A11's `m = 2` σ on p₂ from
0.1786 (8 cells, A) to 0.1500 (9 cells, B). Both readings are "candidate".
These rows are not merged into the F1 files of record, so `analyze
--amendment2` reads A alone.

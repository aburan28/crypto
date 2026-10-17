# XP-20261015-pkm-m3-n21-m4-n20: lead 1 at m = 3, N = 21 and m = 4, N = 20

**Status: blocked, rescoped to `m = 3` (2026-10-16 claim audit).** The
`m = 4, N = 20` cell is **already answered, NEGATIVE**:
[`research/ic_next10_20261014/03-pkm-degree/`](../../research/ic_next10_20261014/03-pkm-degree/README.md)
ran it on 2026-10-14. On both targets msolve 0.6.5 completed three productive
degree-8 F4 rounds (1,241, 2,507 and 3,578 new basis elements) before its 5 GB
cap; because msolve reduces every lower-degree pair first, the degree-7
closure does not contain `1`, so **`D ≥ 8` at `m = 4, N = 20`**, and the
`m = 4` line reads `D = 6, 7, 7, ≥ 8` at `N = 8, 12, 16, 20`. By this card's
own rule that is strong evidence against boundedness at `m = 4`, the only
regime where a tower base could beat rho's exponent
([PKM note](../../research/notes/index-calculus/RESEARCH_PKM_TOWER_ORACLE.md)
§3.7). Caveat: msolve on the raw presentation; `f4_fp_tower` (reduced
presentation) has not run `N = 20`. The card's earlier text calling lead 1
"the only open question" missed this result. The `m = 4` shard is dropped;
the `m = 3, N = 21` cell is still inconclusive there (`D ≥ 7`, stopped by the
wall before any degree-8 round) and is what remains.

Blocked on hardware (a 256 GB host for `m = 3`) or on a disk-backed `B'`
reducer that is not built, and on PR #1893, whose head commit `17c5663` the
card pins. Funding 1,000,000.

## Question

Does the Kummer tower's F4 solving degree `D` stay at 7 at `m = 3, N = 21`,
or does it rise? (The card first also asked about `m = 4, N = 20`; that cell
is answered, `D ≥ 8`; see the status above.)

## Why it matters

Encyclopedia §8 lead 1 asks whether the tower's solving degree is bounded.
At `m = 4` it is answered negatively (above); at `m = 3` it stays open. It
was written as the only open question in the PKM line that decides an *exponent*: a
bounded `D` is what could make a tower base sub-exponential to decompose
against. The other PKM cards (the P-192 auxiliary-curve search, the P-224
depth checks) are worth their cost only if this answer is "bounded".

- `D` rises to 8 at either cell: strong evidence against boundedness at that
  `m`; deprioritise the tower constructions.
- `D = 7` at `m = 3, N = 21`: the plateau from `N = 15` grows to `L = 6`.
  Under rule A1 that is not yet H1a (which needs `L > 10`, `N ≥ 27` at `m = 3`);
  it justifies the larger runs.

## Prerequisites

- A host with 192–256 GB (`m = 3`) or about 512 GB (`m = 4`) of RAM, or a
  disk-backed `B'` streamed from NVMe (not built; would need 130–380 GB of
  disk and one read of `B'` per pass; I/O cost not measured).
- `pkm_tower_pilot` at `17c5663` (PR #1893), or that PR's merge commit.
- `PROTOCOL.md` in this directory, frozen before the card is ready: the
  `--budget`/`--max-nnz` caps (or `--max-dense` and the address-space cap),
  targets per cell, and the checks. `pkm_tower_check verify` can search
  exhaustively at both cells:
  - `m = 3, N = 21` (`t = 7`, `|V| = 128`) is `128³ ≈ 2.1·10⁶` triples per
    row;
  - `m = 4, N = 20` (`t = 5`, `|V| = 32`) is `32⁴ ≈ 1.0·10⁶` quadruples per
    row;
  - it has already checked finished rows at `m = 4`, `N = 16` (round 3: 2
    rows, both agree) and at `m = 3`, `N = 18` (round 7).

  No separate independent check is needed for the verdict. The `m = 3` cell's
  protocol is
  [`research/pkm_lead1_m3n21_20261016/PROTOCOL.md`](../../research/pkm_lead1_m3n21_20261016/PROTOCOL.md),
  frozen; the `m = 4` cell has none yet.
- Resumability: the engine has no checkpoint, so a run that dies loses
  everything and the card is not preemptible. Ready requires either a
  checkpoint or a protocol that accepts that.

## Expected cost (derived in PR #1893, §5; not run)

| cell | growth assumption | width | memory | multiply-adds | time, one thread |
|:--|:--|--:|--:|--:|--:|
| `m = 3, N = 21` | last ratio | ≈ 360,000 | ≈ 146 GB | ≈ 2.5·10¹⁵ | ≈ 86 h |
| `m = 3, N = 21` | geometric mean | ≈ 475,000 | ≈ 133 GB | ≈ 7.4·10¹⁵ | ≈ 138 h |
| `m = 4, N = 20` | last ratio | ≈ 625,000 | ≈ 384 GB | ≈ 8.5·10¹⁵ | ≈ 327 h |
| `m = 4, N = 20` | geometric mean | ≈ 658,000 | ≈ 243 GB | ≈ 2.0·10¹⁶ | ≈ 168 h |

The `m = 4` rows are kept for the record and will not run: that cell is
answered (`D ≥ 8`, above), so its projection, which assumed `D` stays at 7,
was too low. The inputs are the recorded Kummer rows at `p₁ = 15·2^27 + 1`. Each
projection applies one more growth step and assumes `D` does not rise; if it
rises, the projections are too low. The round-6 full-rank exit cuts the time
to about 0.61× and leaves memory unchanged. The kernel is rayon-parallel, so a
many-core node shortens wall time roughly by its core count.

## What "done" means

Both cells ran to a verdict (`D` reached, or a recorded cap), with every row,
the run receipts and the independent check landed in
`research/pkm_lead1_m3n21_m4n20_<date>/` by a PR that updates encyclopedia
§8 and the PKM note. A timeout or an out-of-memory kill is recorded as such
and is not evidence about `D`.

Execution: two cairn lab tasks, one per cell; see [backlog/README.md](../README.md).

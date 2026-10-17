# Generic-DLP correctness fixes: collab 2-cycle escape, Pohlig–Hellman

2026-10-17. Protocol: [`PROTOCOL.md`](PROTOCOL.md), frozen before any count
was taken. Report: [`REPORT.pdf`](REPORT.pdf) (source `REPORT.typ`), with
the figure [`fig_correctness.svg`](fig_correctness.svg). Every number here is
a deterministic count read from [`summary.json`](summary.json), which
`examples/generic_dlp_correctness_summary.rs` writes from `raw/` and
`transcripts/`. Nothing is timed.

These are findings 3–6 of the generic-DLP audit,
[`research/generic_dlp_memory_parallelism_20261017/AUDIT_FINDINGS.md`](../generic_dlp_memory_parallelism_20261017/AUDIT_FINDINGS.md).

## What changed

1. **`pollard_collab` negation-map 2-cycle escape**
   (`src/cryptanalysis/pollard_collab/walk.rs`).
   - The 2-cycle test compared the *unflipped* sum with the *normalised*
     previous point, so it never fired. Walkers caught in a fruitless
     2-cycle ran to the `20·2^dp_bits` step cap and were abandoned.
   - It now compares after the flip and escapes by doubling the smaller-x
     member of the cycle, then normalises again.
   - The escape depends only on the walk's current and previous points, so
     walkers that have merged stay merged.
2. **Pohlig–Hellman digits by baby-step giant-step**
   (`src/cryptanalysis/pohlig_hellman.rs`).
   - Each base-`q` digit was a linear search over `[0, q)`.
   - Now one `m = ⌈√q⌉`-point baby table per prime serves every digit of
     that prime, and a digit costs at most `m` giant steps.
   - Each digit is still the smallest solution below `q`, so every residue
     and answer is unchanged.
3. **Pohlig–Hellman factorisation bounded.** `factorise_up_to(n, bound)`
   stops trial division at the smoothness bound. The old `factorise_small`
   divided by every integer up to the square root of the remaining
   cofactor, about `2^128` divisions on a 256-bit prime order.
4. **Pohlig–Hellman answers checked.** The CRT answer is returned only if
   `[d]G = Q`. A generator of order 1 yields `d = 0` when `Q = O`.

## Arms

- **Before:** `main` at `b1bf5b893`, with this study's driver
  `examples/generic_dlp_correctness_counts.rs` added. The driver uses only
  API present on `main`.
- **After:** the fix commit, first recorded on `b1bf5b893` and then rebased
  onto `main` at `89036527b`.
  - `main` did not change any of the touched sources in between: `walk.rs`,
    `pohlig_hellman.rs`, `ec_vartime_walks_pinned.rs` and the `ecc/`
    stack.
  - The driver re-run on the rebased build reproduces
    `raw/counts_after.jsonl` byte for byte.
- `raw/manifest.json` records the two binaries' SHA-256 values and the
  toolchain, rustc 1.97.0.

## Curves

All four curves are registered (`inputs/curves.json`). PARI/GP 2.15.4
checked each one's group order (`ellcard`) and generator (`ellisoncurve`,
`ellorder`).

| name in code | ICV1 identity | used for |
|:--|:--|:--|
| `demo-small` (`demo-10007`) | `ICV1:fp-10007:-31:10039:6177:unk:unk:r:68b86b706924` | collab |
| `demo-mid` (`demo-99013`) | `ICV1:fp-99013:121:98893:1152:unk:unk:r:95219f0d1db6` | collab |
| `demo-32` | `ICV1:fp-3094378543:-101415:3094479959:3042964895:unk:unk:r:9b187f52394e` | collab |
| `pinned-j0-smooth47` | `ICV1:fp-140737529461027:-2395735:140737531856763:0:unk:unk:r:7a78ff243e18` | Pohlig–Hellman; order `3·19⁴·151·181·13171` |

## Results

### Collab walkers

Target `Q = ⌊2n/3⌋·G`, job seed `0xC011AB`, the default `dp_bits` of
`JobSpec::new`.

| curve | negation map | walkers | DPs before → after | dead trails before → after | DP records verifying |
|:--|:--|--:|:--|:--|:--|
| demo-small | no | 512 | 512 → 512 | 0 → 0 | all |
| demo-small | yes | 512 | 462 → 511 | **50 → 1** | all |
| demo-mid | no | 512 | 511 → 511 | 1 → 1 | all |
| demo-mid | yes | 512 | 406 → 512 | **106 → 0** | all |
| demo-32 | no | 256 | 256 → 256 | 0 → 0 | all |
| demo-32 | yes | 256 | 59 → 241 | **197 → 15** | all |

- **Without the negation map nothing changes**, as the protocol required:
  the counts, steps and DP records are all identical.
- **With it**, dead trails fell from 353 to 16 across the three curves.
  Every DP record passes `DpRecord::verify` in both arms.
- **The 16 that still die** were not classified walker by walker. Longer
  fruitless cycles (4, 6, … steps), which this fix does not handle, are the
  expected cause, and the step cap still bounds them.

### Pohlig–Hellman

24 targets (xorshift64 seeded `0x9E3779B97F4A7C15`, reduced mod `n`) × 2
bounds (`2^14`, `2^17`), 48 solves in all.

- All 48 are correct before and after, and every recovered logarithm is
  identical.
- `total_steps` went from 284,364 to 10,682 (**26.62×**). Per solve the
  ratio is 1.85–46.53×, and it is above 1 on every solve.
- The two bounds give identical counts, because every prime factor of `n`
  is below `2^14`.
- The low end of the range is where a target's digits are small: the linear
  search stopped early there, while baby-step giant-step always pays for
  its `⌈√q⌉` table.

## Pinned walk test re-pin (`tests/ec_vartime_walks_pinned.rs`)

The collab fix changes walk transcripts by design, and the Pohlig–Hellman
change changes `total_steps`, so both pinned digests change.
`transcripts/` holds the full transcripts the test hashes, before and after.
The "before" transcripts came from `main`'s sources with the expected
digests zeroed, and they reproduced the pinned digests exactly. The old
digests are kept in comments beside the new ones.

- **Collab.** The `neg=false` sections are byte-identical. The `neg=true`
  sections went from 30 dead trails to 16.
- **Pohlig–Hellman and Cheon.** Apart from step counts, the two transcripts
  are identical: every answer, factorisation, residue and wrong-order
  outcome, and every Cheon line. Over 72 Pohlig–Hellman lines,
  `total_steps` went from 218,517 to 8,648.

## Class and scope

- **Class.** The collab change is a correctness fix. The Pohlig–Hellman
  change is an algorithmic change to its counted cost, from `O(e·q)` to
  `O((e+1)·√q)` per prime power `q^e`.
  - Neither touches `S`, the index-calculus engines, or any IC/rho
    comparison.
  - No scoreboard, leaderboard or progress figure changes. The curve
    registry and its derived files were regenerated only to add the four
    study curves.
- **What the counts cover.** They cover these curves and targets only. How
  often fruitless 2-cycles arise elsewhere depends on the branch count and
  `dp_bits`.
- **What is not measured.**
  - Pohlig–Hellman wall time; the counts are the metric.
  - An external comparison (§15). Reducing to prime-order subgroups and
    solving each with a square-root method is the standard Pohlig–Hellman,
    as in PARI/GP `elllog` and Sage's generic `discrete_log`. The change
    brings this module up to that method; it does not claim to go past it.

## Verification

- `cargo test --release --lib`: 4232 passed, 0 failed, 117 ignored, on the
  branch before the rebase onto `89036527b`.
- New tests:
  - `fruitless_two_cycles_are_escaped_deterministically` (`walk.rs`);
  - `large_prime_order_stops_at_the_bound` and
    `answers_are_verified_in_the_group` (`pohlig_hellman.rs`).
- The existing differential test now compares answers only and bounds steps
  by `m·(e+1)`.
- **Known failure on `main`, not this change.**
  `tests/ecdlp_nist.rs::list_covers_the_fifteen_nist_curves_and_finds_no_structural_weakness`
  fails on `main` too: it counts 24 curves where it expects 15, since
  `46481bca5` added the Certicom curves. CI only runs `cargo check` on test
  targets, so it does not catch it.

## Files

- `PROTOCOL.md` — the frozen protocol.
- `inputs/curves.json` — the study curves.
- `raw/counts_{before,after}.jsonl` and `raw/manifest.json` — the counts and
  the run manifest.
- `transcripts/` — the pinned test's transcripts.
- `summary.json` — native summary.
- `scripts/plot.py` — layout only.
- `fig_correctness.{svg,png}`, `REPORT.typ`, `REPORT.pdf`.

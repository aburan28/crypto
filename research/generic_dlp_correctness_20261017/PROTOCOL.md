# Protocol: generic-DLP correctness fixes (collab 2-cycle escape, Pohlig–Hellman)

Frozen 2026-10-17, before any count was taken. Study:
`research/generic_dlp_correctness_20261017/`. Findings 3–6 of
`research/generic_dlp_memory_parallelism_20261017/AUDIT_FINDINGS.md`.

## Question

The audit found four correctness defects. This study asks whether the fixes
behave as intended, in deterministic counts. Nothing here is timed.

1. **`pollard_collab` negation-map 2-cycle escape.** The 2-cycle check
   compared the unflipped sum with the normalised previous point, so it never
   fired.
   - The fix compares after the flip and escapes by doubling the smaller-x
     member.
   - Hypothesis: with the negation map, fewer walkers end as dead trails at
     the step cap, more reach a distinguished point, and every DP record
     still verifies.
   - Without the negation map nothing changes: the code path is not taken.
2. **Pohlig–Hellman digits.** They were linear search, `O(q)` per digit.
   They are now baby-step giant-step, `⌈√q⌉` table points per prime and at
   most `⌈√q⌉` giant points per digit.
   - Hypothesis: the same recovered logarithm on every target, with
     `total_steps` lower.
3. **Pohlig–Hellman factorisation.** It is now bounded by the smoothness
   bound; previously trial division ran to the square root of the cofactor.
4. **Pohlig–Hellman answer.** It is now checked in the group (`[d]G = Q`),
   and order 1 yields `d = 0`.

## Arms

- *Before*: `main` at `b1bf5b893`.
- *After*: the branch `claude/generic-dlp-correctness-fixes`.
- The same driver, `examples/generic_dlp_correctness_counts.rs`, uses only
  API present on `main`, so it builds on both.

## Frozen inputs

- **collab.** `demo_curve("demo-small")`, `demo-mid` and `demo-32`.
  - Target `Q = ⌊2n/3⌋·G`; job seed `0xC011AB`.
  - The default `dp_bits` of `JobSpec::new`.
  - Walkers `0..512` (`0..256` on `demo-32`), with and without the negation
    map.
- **Pohlig–Hellman.** The 47-bit `y² = x³ + 4` curve of
  `tests/ec_vartime_walks_pinned.rs`, of order `3·19⁴·151·181·13171`.
  - 24 targets from xorshift64 seeded `0x9E3779B97F4A7C15`, reduced mod `n`.
  - Bounds `2^14` and `2^17`.

## Measures

- **collab.** Per curve and mode: DPs reached, dead trails, steps of each,
  and how many DP records pass `DpRecord::verify`.
- **Pohlig–Hellman.** Per target and bound: the recovered logarithm, whether
  it equals the planted one, and `total_steps`.

## Success, rejection

- **Reject** the fix if any recovered logarithm differs from the planted
  one, any DP record fails verification, or anything changes in the
  no-negation collab arm.
- **Success:**
  - (1) negation-map dead trails fall;
  - (2) PH `total_steps` falls on every target, with the same answer;
  - (3) and (4) hold by their tests.

## What this cannot establish

- No wall time and no `S`. The collab change alters walk transcripts by
  design, so the pinned digests in `tests/ec_vartime_walks_pinned.rs` are
  re-pinned deliberately, with the before/after counts as the record.
- No claim about how often 2-cycles arise in another configuration.

# XP-20261015-kummer-p224-depth: does the Kummer solving degree transfer to the real field?

**Status: blocked** on the generic-field F4 lane (`claude/f4-generic-field`,
in progress, not pushed). Funding 400,000, second of the four seed cards.

## Question

Every recorded PKM solving degree was measured over a prime below `2^32`
(mostly `p₁ = 15·2^27 + 1`). Is `D` at fixed `(m, N)` the same at a 45-bit
prime with `v₂(p − 1) = 40` (`p = 27·2^40 + 1 = 29686813949953`, prime by GNU
`factor`) and at the NIST P-224 prime `2^224 − 2^96 + 1` (`v₂(p − 1) = 96`)?

## Why it matters

- `D` equal to the `p₁` value in every cell: lead 1's answer from the cheap
  `p₁` runs speaks for P-224, and the small-prime ladder stays the instrument.
- `D` different in any cell: every `p₁` extrapolation in the encyclopedia is
  suspect, and lead 1 must be rerun at the real field.
- Towers deeper than 27 levels exist only at primes like these, so this card
  is also where a deeper Kummer tower can first be built.

## Prerequisites

- The tower engine over a 45-bit and a 224-bit prime (the generic-field lane).
  Its flags may differ from today's `pkm_tower_pilot`; `PROTOCOL.md` fixes the
  final command.
- `PROTOCOL.md` here, frozen before ready: the cells (the card repeats
  `m = 2, N ≤ 20`; `m = 3, N ≤ 15`; `m = 4, N ≤ 12`), targets per cell, and the
  comparison with the recorded `p₁` rows.

## Expected cost

`cost = (recorded p₁ cost of the cell) × F`, with `F` the multi-limb slowdown
of the generic-field engine, which that lane measures. The cells are small at
`p₁` (`m = 3, N = 15`: 1.1 GB, 6.0 minutes; the others less), so with `F` in
the tens the card is hours to a few days of one core. Memory grows with the
limb count: about ×2 at 45 bits and ×7 at 224 bits for 32-bit limbs.

## What "done" means

Every cell at both primes reached `D` or a recorded cap, compared row by row
with the `p₁` values, landed in `research/pkm_kummer_wide_field_depth_<date>/`
by a PR that updates encyclopedia §0.1. Execution: two cairn lab tasks, one
per prime; see [backlog/README.md](../README.md).

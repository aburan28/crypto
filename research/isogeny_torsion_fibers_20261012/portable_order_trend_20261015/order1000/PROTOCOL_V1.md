# Extended different-j prime-to-characteristic order test

**Frozen before running the screen or any child solve.** This is a
bounded extension of the already checked order-229/421/433/457 panel,
using the same source family `y²=x³−21x+37`, rational order-three
kernel `T=(3,1)`, explicit Vélu map and deterministic Rust arithmetic.
It tests whether the counted full-image/rho trend persists at an order
above 1,000. It does not infer an asymptotic exponent from five cells.

## Deterministic selector

Run native `examples/different_j_extended_screen.rs screen-larger` once.
Source SHA-256:
`272aeb4ac002ef8c623564bf8eebc8b5b29e97592afe23c31fc433d769d07b3f`;
release binary SHA-256:
`ea6e3166529f22d23b5373acb55f8b8fc6827b5d3ef928a8ccab57b3d32bc434`.
Enumerate every prime `p` in ascending order from 619 through 1601 with
`p≡1 (mod 3)`. Require nonsingular source/destination, distinct `j`,
matching base-field order, the kernel identity and 64 seeded independent
map-law checks. Factor `#E(F_(p²))`; accept a prime divisor `r≥1009` of
multiplicity one with `r≠p`. The first accepted row in ascending `p`
is the **single** follow-on cell, regardless of its apparent performance.
Preserve all rejects and the summary. Whole screen: 600 seconds and 512
MiB resident ceiling. If no row qualifies or any check fails, retain
the failed run and stop before target generation.

## Frozen follow-on comparison

For the selected `(p,r,h)`, use the existing native
`different_j_trend_fiber` exhaustive source/image/preimage census and
independent `different_j_trend_fiber_check`, including the entire
`r²` subgroup pair map-law check. Require every rational image fiber
count and exact order `r`. Freeze 64 distinct public source targets
with `different_j_trend_batch freeze`; keep scalar witnesses checker-only.
Commit source, binary, public-input, checker-input and fiber receipt
hashes in a run spec **before** child measurement. No replacement
target or changed cap after inspection.

Run the same seven arms, three rotated rounds and prefixes
`N=1,4,16,64` as the V2 panel: source projected, isogeny image,
destination direct, full rational image, source BSGS, source negation
rho and rho A/A. Each arm/round/prefix is a fresh cold process,
so 84 children are expected. Keep eight signed columns, relation cap
`8r` trials per target, rho at most 32 restarts of `20r` steps, and
BSGS baby size `ceil(sqrt(r))`. Each child has a 300-second wall
limit and 512 MiB peak resident ceiling. Preserve capped, timed-out,
or failed children; do not impute missing scalars.

Use `different_j_trend_batch_check` to independently rebuild all
field/curve/map/row arithmetic, every relation witness and every
`[k]G=Q` and image identity on all returned scalars. Report the exact
number checked; a failed replay blocks a complete-solve result.
Compare cold setup through verification group counts and the separate
Fp multiplication/inversion axes against direct source BSGS and rho
on identical target prefixes. Report one `S_group=(adds+doubles)/
(N sqrt(r))` table, with rho A/A, source and destination relation arms,
all costs and memory. Do not add overlapping field and group counts.
Complete portable S remains unset until `../PRICE_PROTOCOL_V2.md` is
satisfied. L2 wall timing and independent-host performance are distinct
gates.

## Decision and provenance

The bounded trend is supported if the selected cell completes all 84
children, independent replay agrees, and every frozen count is
reported. It is contradicted or unresolved on its measured axis if
full-image arithmetic falls below rho, a cap is reached, or a child
fails. A selected-cell outcome does not establish a general isogeny
advantage or a fitted exponent. The native screen's exact source and
binary hashes, environment, command, raw rows, and exit status belong
beside this protocol. Visual run/results PDFs and benchmark mirror
follow the repository evidence rules.

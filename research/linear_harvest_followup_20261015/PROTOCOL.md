# Linear-harvest follow-ups (DRAFT, not frozen)

Status: **planned, no data collected.** Tracking PR. Follows
`research/ic_improve_20261015/05-linear-harvest/` (round 2 item 05) and its open follow-ups.

## Questions
1. **End to end at (23,12).** The bench gain (1.166x) did not carry through because the
   production path stops at the first accepted root. Does the gated root harvest change
   production solve counts or counted cost at (23,12)?
2. **Gate on variable count or tree size.** The gpv < 0.9607 gate skips the largest gain
   (ungated root harvest 1.28x at (29,15) m2). Is there a gate on variable count or tree size
   that fires there without regressing any training cell?

## Design (to freeze before any data)
- Incumbent: harvest off. Candidate: root-only harvest with the new gate. Same targets, same
  counted unit, rotated arm order, every root set verified.
- Choose the gate on training targets only (t = 0..7 style split as in item 05), freeze, then
  test on held-out targets; admissible = no training cell with I/G CI upper bound below 1.
- Level of the change stated per the ecbench-bounds vocabulary (ADVANCES / TRADES / MATCHES /
  REGRESSES; constant vs exponent).
- Counted operations or thread CPU only, not wall-clock on the shared host.

## Out of scope
Changing the default gate in `koblitz_groebner.rs` before a held-out ADVANCES verdict.

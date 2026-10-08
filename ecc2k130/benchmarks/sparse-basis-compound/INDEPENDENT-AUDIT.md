# Independent audit of the sparse-basis compound screen

Audited commit: `cfa375b456fa43dc566e1aef8a6d63bb65b81bd2`.

Decision: **PASS with the limitations below.**  A fresh native compilation
reproduced `result.json` byte for byte at SHA-256
`9608c61a06179bcc6be63e0892c061a84c9217f872f711b4bb16ebe29f5ff66b`.
Every entry in `files.sha256` verified, the worktree stayed clean, and no GPU
was launched.

## Recomputed findings

- Both pentanomials satisfy the prime-degree Rabin conditions.
- `ord_263(2) = 131` and `2^131 mod 263 = 1`; both reported beta roots
  evaluate to zero in their candidate fields, and both power matrices have
  rank 131.
- Per candidate, the native audit covers all 17,161 field-basis product pairs,
  all 261 raw reducer directions, all 131 conversion basis vectors, and all
  1,179 normal/`L_j` basis-map directions.  Reduction and the maps are linear,
  and multiplication is bilinear, so these basis panels prove the corresponding
  full-space identities.  The deterministic dense panels are extra
  implementation guards.
- Every `L_j`, `j=3..10`, has rank 130, as required by `gcd(j,131)=1`.
- The sparse two-fold reducer recounts to 40 shifts, 16 ORs, 24 XORs, and two
  ANDs: 82 source word operations.  The current direct reducer recounts to
  86 shifts, 39 ORs, 53 XORs, and two ANDs: 180 operations.
- The corrected batch ledger is 85 products per 16 updates.  The 31 reverse
  products already include the 16 `inv*W_i` lambda products.  With
  `PACKED_INV_POLY=0`, the inverse's eight products use
  `fromPolynomialProduct131` and do not call the direct beta reducer.
- The resulting counts are 5.8125 sparse reductions and 38.125 `CLMAD`
  lane-instructions/update.  The direct-reducer subledger is 1,046.25 current
  versus 476.625 sparse source operations/update.  The conservatively charged
  reducer-plus-map saving is 809 operations/update.
- Eleven maps occupy 77,440 bytes, or 78,464 bytes with the reservation.  At
  batch 16 they cost 275 LDS instructions, 687.5 scalar-word equivalents,
  984.375 exact evaluator ALU operations, or 990.625 under the frozen
  conservative charge.
- The one-pipe arithmetic recomputes to 15.283374545 B/s for the transferred
  random-`LDS.U8` model, 26.579781818 B/s at the 16-lane sensitivity,
  29.514458044 B/s for charged table ALU, and 23.965377049 B/s for `CLMAD`.
  The frozen 5% admission threshold is 16.208510850 B/s.

The 15.436677 B/s reference was independently read from commit
`e0b0858b179fe28f753353e7a7313cfe23d71884`.  The headline result and audit
match the pinned SHA-256 values in the protocol.

## Limits on the decision

The 15.283375 B/s lookup figure transfers a measured random-`LDS.U8` rate to
an aligned eight-entry `LDS.128` plus `LDS.32` multicast pattern.  It is not a
hardware upper bound for that exact pattern.  The 16-lane sensitivity passes
the admission threshold, so the exact-layout microprobe can reopen the route.

The register counts are source-level overlap estimates.  Only a CUDA compile
can establish allocation, spills, and whether a 512-thread block launches.
The source-operation ledger is not a SASS prediction and omits unchanged
inverse transforms, address generation, state traffic, loop/control work, and
initialization.  Those omissions favor the candidate in the rejection model.

Accordingly, `NEGATIVE_STATIC_EVIDENCE_UNDER_FROZEN_TRANSFER_MODEL` and “do
not build the whole-walk candidate under this protocol” are supported.  A
broader claim that sparse-basis tables are impossible on this GPU is not.

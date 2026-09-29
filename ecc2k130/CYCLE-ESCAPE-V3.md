# Table-walk cycle escape v3

Rule v3 fixes the entry-dependent exits found in the September 21 review.
The intervening v2 change already added both four-term Frobenius relations
and pairwise exclusions covering cycles through length eight. That predicate
is retained as a cheap **hint**, not permission to change a step.

On a hint, host replay and packed selection run the same bounded procedure:

1. Follow at most eight **raw** table steps, without history-based changes.
2. Require an exact return to the starting point. An open path is not a cycle.
   If the probe sees a distinguished point, keep the raw step so it is reported.
3. Recompute each cycle vertex's eligibility from its four cyclic predecessor
   tags. Choose the least eligible `(normal-basis weight, canonical x)` key.
4. Advance the branch once only at that anchor; elsewhere keep the raw step.

Thus the incoming history may delay discovery, but cannot choose the exit.
For a detected cycle, a cold arrival gets a full four-tag history within four
steps and reaches the common anchor within another lap. Frobenius and negation
preserve the ordering and the eligibility test, so related cycles use related
exit edges. This is eventual coalescence with a bounded delay, not equality of
the very next step for every pair of different histories. A finite trail guard
may still truncate a delayed arrival; no claim of zero campaign collision loss
is made.

The regression includes the actual F131 raw two-cycle containing `[1184]P`,
whose tags are `0x03e0` and `0x13e0`. Both entry phases and their Frobenius and
negation images must leave through one orbit. The packed selector is compared
with the reference on that cold path, including the hybrid table layout.
Synthetic tests cover all entry phases of pairwise cycles through eight steps,
both four-term tau families, false hints, and cycles containing a DP.

## Limits and cost

The probe detects exact point cycles through eight steps when the tag predicate
provides a hint. It does not certify the absence of longer cycles, short orbit
cycles with a longer point period, residual six-term relations without a hint,
or cycles involving the modified escape edges. The standalone default is now
`--max-iters 4294967296`; an explicit `--max-iters 0` disables it. Checks occur
on the existing 4096-step boundaries. Campaign users should size the limit for
their DP distribution; it is not part of the iteration's algebra.

The cold path performs extra additions/inversions and possibly canonicalization.
Its GPU register, occupancy, throughput, and end-to-end cost have **not** been
measured. Existing v2 throughput and walk-constant figures do not transfer to v3.
This repair makes no speedup claim and does not change the default sigma walk.
The hashed controls in `walkconstant.cpp` explicitly retain v2; native rows use
v3 and identify it in their output. Their old merge-loss model is not a v3 model.

## Compatibility

V3 changes the seed-to-endpoint map. Start a new table corpus and retain the old
client for old data. Table DP files now begin with `ECC2KDT3`, version 3 and a
32-byte record size in the 16-byte header. Sigma clients reject that header;
v3 table clients reject old headerless or sigma witness data. New table
checkpoints use reference version 34 and packed version 35. The AWS campaign
contract includes the new rule identity, and spool/merge/ingest understand the
new framing. No live campaign is migrated or restarted by this change.

`make test-cycle-rule test-cycle-escape` runs the hint and escape regressions.
The production test separately checks replayed collision recovery, corpus
isolation, guard behavior, and checkpoint compatibility. See
`research/ecc2k_cycle_escape_20260926/` for the frozen protocol and local receipts.

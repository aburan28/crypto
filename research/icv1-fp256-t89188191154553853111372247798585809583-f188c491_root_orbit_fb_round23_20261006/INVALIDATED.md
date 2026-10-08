# Invalidated round-23 execution

The first release execution on 2026-10-06 stopped during the first arm's
independent point replay with 540 mismatches among 4,096 rows.  It emitted no
result artifact and receives no correctness, cost, or screening credit.

Cause: the independent BigUint path emitted minimal-length big-endian bytes,
while the canonical P-256 rows use fixed-width 32-byte field encodings.  Values
with leading zero bytes therefore compared unequal even though the field
elements agreed.

Correction: the independent path now rejects values longer than 32 bytes and
left-pads shorter values to exactly 32 bytes.  The canonical execution replayed
4,096 sorted columns in each of 16 arms (65,536 total) with zero mismatches.
The stopped execution produced no JSON to preserve or cite.

# Deviations from PROTOCOL_EDWARDS.md (recorded before the affected arms report)

1. **Launch (17:22–17:23Z).** An earlier launch attempt that had seemed not to
   start was in fact running (from 16:47Z). A second launch overlapped it, and
   both appended to `dlp_p31_s1.jsonl`. Both were killed, `data/edwards/` was
   deleted, and one orchestration was started at 17:23Z on commit
   `3a8599b3f`. Nothing from the overlapped launches is kept.

2. **Order of the DLP arms (20:16Z).** `dlp_p41_s1` drew a curve with a
   six-point base (`C(6, 4) = 15`): the floor is about `2.4·10⁵` trials per
   relation, about `1.7·10⁶` trials for 7 relations, roughly 13 h, so it will
   probably hit its 10-hour cap.
   - It is left running to its cap. Its process continues; only the
     orchestration that would have followed it was detached. Because of
     that, a timeout of this arm is not written to `exits.txt`; its record
     then simply ends without an `ic` row.
   - The remaining arms (`p = 61` and `71`, seed 1; `p = 31` and `41`,
     seed 2) run in a second sequence alongside it from 20:16Z, each with
     the same 10-hour cap.
   - Every requested arm still runs; only their order and concurrency
     changed. Counts are unaffected; wall times overlap.

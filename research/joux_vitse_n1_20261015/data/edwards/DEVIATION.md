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

3. **Container restart (about 02:05Z, 2026-10-16).** The restart killed
   `dlp_p41_s1` (started 18:15Z, about 7.8 h, 1 relation by trial 165,664)
   and `dlp_p71_s1` (started 00:02Z, about 2 h, 2 relations by trial
   155,645).
   - Their partial records are kept in `killed_20261016T0205Z/` and are not
     used in any result.
   - The arms still owed by the protocol were relaunched at 02:14Z in one
     sequence, each with its 10-hour cap: `p = 31` seed 2, `p = 41` seed 2,
     `p = 71` seed 1, `p = 41` seed 1.
   - The two seed-1 arms restart from trial 0 and are deterministic, so they
     replay the killed runs' relations first.

4. **Resumable DLP driver (protocol addendum, 2026-10-16 ~03:30Z, before
   the arms it affects start).** After two container restarts, `dlp` mode
   writes a checkpoint after every batch (`<record>.ckpt`), and `--resume`
   continues from it.
   - A resumed run rebuilds the deterministic set-up, skips the rows it
     already wrote, and fast-forwards the target generator past the
     completed batches. It therefore draws the same targets in the same
     order as an uninterrupted run.
   - The record's counts, relations and answer are those of an uninterrupted
     run; only wall time (summed over segments) and a `resume` row differ.
   - Without a checkpoint the run starts fresh and truncates its record
     itself.
   - The trial, the trace, the linear algebra, rho and the accounting are
     unchanged. Equivalence was checked before use: an uninterrupted run and
     a split-and-resumed run of the same arm, compared row for row (see the
     run report).
   - The arm running at the time (`dlp_p31_s2`) continues on the old binary;
     later arms run the new one.

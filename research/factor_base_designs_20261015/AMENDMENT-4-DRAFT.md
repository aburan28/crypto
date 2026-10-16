# Amendment 4 — DRAFT, not registered

**Status:** proposal from the experiment lane, 2026-10-16. It binds nothing
until it is appended, dated, under `## Amendments` in
[`PREREGISTRATION.md`](PREREGISTRATION.md), before any G1 row exists. Any
wording may change before then. The text between the rules below is the
proposed amendment.

---

### Amendment 4, 2026-10-16: the `m = 4` series to three cells; distinct targets

Registered after F1–F4 (`runs/analysis-amendment2.txt`) and before any G1
row exists. It adds one experiment (G1) and one construction rule (G2). It
changes no committed row, verdict, threshold or decision rule. The
instrument, the seed rule (Amendments 1 and 2), the completed-cell rule, the
σ fit and the secondary width slope are unchanged, except where stated.

**Why.**
- F4 left both `m = 4` series inconclusive with 2 completed cells each.
- At `t = 4`, every system stopped at the 300 s budget. A0 stopped at
  partial D = 7 with 102 516 columns, and A11 at partial D = 6 with 88 740
  columns.
- §1's threshold was derived at `m = 4`, so `m = 4` is where a σ verdict
  means most.
- F3 showed that a 3600 s budget can complete a cell that 300 s cannot. A11
  at `m = 3`, `t = 5` completed in about 570 s.
- F3 also showed the limit that may bind next. A0 at `m = 3`, `t = 5`
  exhausted `ulimit -v 13000000` after about 1545 s.

**G1. `m = 4` from `t = 4`, p₁, arms A0 and A11.**
- *Budget:* F4 budget 3600 s and hard kill 3900 s, as in F3. The memory limit
  (`ulimit -v 13000000`) is unchanged. A system that hits it is recorded with
  its reason.
- *Seeds:* as in Amendment 2. The `t = 4` systems are therefore exactly the
  F4 systems that timed out, now measured with the longer budget.
- *Stop rule:* `t` rises from 4 until the stop rule (Amendment 1 item 5),
  per arm.
- *Files:* `runs/G1-{arm}-m4.jsonl`.
- *Series:* each arm's `m = 4` series is its F4 rows at `t = 2, 3` plus its
  G1 rows at `t ≥ 4`. F4's `t = 4` timeout rows stay committed, and the G1
  rows supersede them in the fit. A cell enters the fit only if both of its
  systems complete.
- *Worst-case cost:* 4 systems × 3900 s ≈ 4.3 h at `t = 4`, plus `t = 5` for
  an arm whose `t = 4` cell completes.

**G2. Distinct targets within a cell (construction).**
- When a planted target's `x_R` equals the `x_R` of an earlier target in the
  same cell, it is rejected. The next draw comes from the same target
  stream, as a failed plant is redrawn.
- *Why:* for the tower arms the system depends on the target only through
  `x_R`, so such a cell holds one system twice. Ten committed cells have it:
  eight on p₁ (A0 m = 2 t = 2, 3; N m = 2 t = 2; A1z m = 2 t = 2, 3 and
  m = 3 t = 2; A1g m = 2 t = 4; A4 m = 3 t = 2), and A11 m = 2, 3 at t = 2 on
  p₂. N11 m = 2 t = 2 (F2) also has it.
- *Scope:* the rule changes a draw only where the old rule collides. Every
  committed row without a collision is reproduced unchanged. It applies to
  G1 and to every later run. Committed collided cells are not re-run;
  `analyze` marks them as single-system cells.
- *Implementation:* target `k` re-plants targets `0 … k−1` from their own
  streams to learn their `x_R`, which is cheap. A regression check
  (`regress`) on committed cells without a collision must still pass before
  any G1 row.

**Not proposed: more targets per cell.**
- Every completed cell so far has equal D across its two systems: 97 cells
  over both primes and all arms, of which 76 have distinct `x_R`.
- 96 of those 97 also have equal `max_cols`. The exception is A1g at
  `m = 3`, `t = 3`.
- At these sizes, instance-to-instance variation in D is not observed, so
  further targets would add cost and little information. The scarce
  quantity is completed cells at larger `t`.

**Decision rules, fixed now.**
- *G1:* the σ verdicts of §1 per arm at `m = 4`, over all of the arm's
  completed `m = 4` cells (F4 and G1).
- *G1, below 3 cells:* if an arm still has fewer than 3 completed cells, its
  verdict stays inconclusive and the reason is recorded (budget or memory).
  No further budget or memory change is made without a new amendment.
- *G1, secondary (descriptive):* Δlog₂ `max_cols`/Δb per arm, and the ratio
  of A11's cell `max_cols` to A0's at each shared completed `t`.
- *G2:* none. It is a construction rule. The rows record `distinct_targets:
  true`.

**Predictions.**
- G1: A0's σ ≥ 0.5 at `m = 4`. This is F4's registered prediction, still
  undecided.
- G1: no prediction for A11's σ. Its `max_cols` stays below A0's at `t = 4`,
  but the ratio is closer to 1 than at `t = 3`.
- G1: A0's `t = 4` systems complete within 3600 s and 13 GB. This one is
  uncertain: F3's A0 cell exhausted memory at a smaller width.

---

**Open for whoever registers this.**
1. **G1's memory limit.** Keep 13 GB, as drafted, or raise it to fit the
   host. A higher limit is a different instrument limit and would be
   recorded in every row. D does not depend on the limit once a system
   completes.
2. **G2's rejection key.** For A11 the system depends on `y(−R)` as well as
   `x_R`. Two A11 targets with equal `x_R` and opposite `y(−R)` are distinct
   systems, yet are rejected by the rule as drafted. Keying on the pair
   (`x_R`, `y(−R)`) for A11 and N11 is the alternative.
3. **The A11 width prediction** is the drafter's reading of the F4 trend
   (A11 at 2159 vs A0 at 21 619 columns at `t = 3`; 88 740 vs 102 516 at the
   `t = 4` stops). Keep it, sharpen it, or drop it.

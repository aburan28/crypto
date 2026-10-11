# Run report: E4, two-torsion quotients of the S₃ system (2026-10-11 17:14 UTC)

Descriptive record of one run (AGENTS.md "Research searches must leave
visual reports").  Interpretation is in the round's
[README](../../README.md).  PDF: [`report.pdf`](report.pdf).

| item | value |
|:--|:--|
| lead | E4 of [`PROTOCOL.md`](../../PROTOCOL.md) |
| command | `isogeny_gap_two_torsion_f4 --pairs research/isogeny_conductor_gap_ic_20261010/pairs/prime_pairs.json --out research/isogeny_conductor_gap_round3_20261011/e4/records.jsonl --sizes 32,64 --budget-ms 20000` |
| runner | `examples/isogeny_gap_two_torsion_f4.rs` at commit `fb8185437`, plus the orbit fix committed with the records |
| solver | `cryptanalysis::f4_fp::f4`, grevlex, degree bound `2·deg + 8`, 20 s budget per system |
| host | Claude Code cloud container (4 vCPU); not an isolated ecbench host, so wall times are observations only |
| started / finished (UTC) | `e4/started_utc` 2026-10-11T17:14Z / `e4/finished_utc` 2026-10-11T17:15Z |
| raw records | [`../../e4/records.jsonl`](../../e4/records.jsonl), run log [`../../e4/run.log`](../../e4/run.log) |
| tables | [`../../TABLES.md`](../../TABLES.md), [`../../tables.json`](../../tables.json), from `examples/isogeny_gap_round3_tables.rs` |

**Workloads.**
- 24 round-1 prime pairs with even cofactor (22 and 26 bits): crater, floor
  and isomorphic floor model, 72 curves.
- Base sizes `B` = 32 and 64 abscissae.
- Subgroups: `1` (direct), `C2`, `C4`, and `E[2]` where rational.
- Three targets per (curve, `B`, `H`): two pseudo-random, one planted sum.

**Requested and completed.**  1,128 systems were requested and 1,128
completed: 752 random-target and 376 planted-target.
- No timeout.
- No system where F4's consistency verdict disagrees with direct
  enumeration of base pairs.
- No planted solution fails to annihilate its basis.
- The Vélu map passed its 32 random-sum checks on every (curve, subgroup).

A first launch looped while building a `C4`-stable base on the first curve.
An orbit repeated an abscissa, so the base could not reach exactly `B`.  It
was stopped, the runner was fixed to skip such orbits, and the run was
restarted from an empty record file.

**Rational 2-torsion found** (recorded per curve, `two_rank`):
- Rank 2 on the crater and rank 1 on the floor and model:
  - all 12 D = −7 pairs;
  - the D = −11 and D = −3 pairs with ℓ = 2.
- Rank 1 at all three roles on the D = −8, −4 and the other D = −3 pairs.

**Measured values** (exact counts, deterministic; the range is over every
curve, role and target).  F4's widest matrix and multiplications depend
only on the base-polynomial degree:

| base-polynomial degree | solving degree | widest matrix (cols) | row-reduction mults, random target | planted target |
|--:|--:|--:|--:|--:|
| 8 | 10 | 45 | 7,408–7,414 | 13,903–13,929 |
| 16 | 18 | 93 | 34,152–34,158 | 81,128–81,178 |
| 32 | 34 | 189 | 145,630 | 483,552–483,561 |
| 64 | 66 | 381 | 600,510 | 3,123,348–3,123,354 |

Derived (shown in `TABLES.md`): direct ÷ quotient at matched `B`.
- Solving degree 1.89 to 3.67.
- Multiplications 4.12 to 4.26 for |H| = 2, and 17.6 to 19.6 for |H| = 4.

On every rank-2 crater, `E[2]` and the floor's and model's `C4` record the
same degree and widest matrix.  Their multiplications agree within 6 of
34,158 at `B` = 64.

**Graphs.**
- [`e4_ops_vs_degree.svg`](e4_ops_vs_degree.svg): every system's
  multiplications against its base-polynomial degree.
- [`e4_gap_b64.svg`](e4_gap_b64.svg): crater `E[2]` against floor and model
  `C4` at `B` = 64.
- Source: [`plot.py`](plot.py), which plots recorded values only.

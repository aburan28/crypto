# ECDLP methods priced in one unit: bit operations

**Report:** [`REPORT.pdf`](REPORT.pdf) (source [`REPORT.typ`](REPORT.typ)).
**Protocol:** [`PROTOCOL.md`](PROTOCOL.md), frozen before any session, with
its deviations.
**Unit:** [`docs/units/BITOPS.md`](../../docs/units/BITOPS.md).

Every ECDLP method in the repository is now charged in `bitop/v1`: two-input
gates per primitive, under frozen formulas. Every solver and pipeline step
counts its work, and an instruction audit checks that the counters see it.
What the study measured:

- **The instrument works.** On a prime ladder of 11 curves (16 to 26 bits,
  528 verified runs), rho, BSGS and kangaroo fit exponents in additions of
  0.433 [0.320, 0.538], 0.474 [0.436, 0.508] and 0.472 [0.370, 0.567]. All
  three intervals contain the theoretical ½.
- **Index calculus with a fixed base scales worse than rho at these sizes.**
  - MITM with a 32-point base: 0.766 [0.600, 0.909].
  - Root-finding with m = 3: 0.925 [0.873, 0.973].

  At 26 bits they cost 30× and 112× the generic floor; rho with negation
  costs 1.42×.
- **The old unit understated algebraic index calculus by up to three orders
  of magnitude.** Each pair is the same execution's committed ratio to the
  floor, then its bit-op ratio:
  - F4 descent: 8.0 → 7,147.
  - SAT: 15.5 → 16,552.
  - Buchberger: 53 → 67,837.
  - inherited-F4: 66 → 144,448.

  Table oracles barely move (within 15%). Signed-Frobenius rho doubles, from
  2.40 to 5.38, because its canonicalisations had been free.
- **Amortization gives no window.** Against batch rho on the prime ladder, no
  index-calculus arm has a window of target counts in which it wins, at any
  size. Every further target costs it 3.5–311× what rho pays.
- **Exponents carry the digits their intervals allow.** About ±0.1 from 88
  runs. A third decimal place would take about 10⁴× the runs.

Wall time is not used. Sessions ran at isolation L0, because bit-operation
counts do not depend on the host. No security claim is made about any
deployed curve.

| directory | what |
|---|---|
| `specs/` | the frozen specs |
| `sessions/` | the sessions, each `ecbench verify --replay 2` clean |
| `bounds/` | `ecbench bound fit --unit bitop/v1` records |
| `audit/` | `ecbench audit-ir` panel, final binary |
| `amortize/` | `ecbench amortize` outputs |
| `archive/` | `ecbench work --summary` over every committed session |
| `figures/` | SVG and PNG figures and `plots.py`, which reads only recorded values |

## Canonical graphs

These canonical graphs were checked and left unchanged:

- `docs/index-calculus-scoreboard.html`
- `docs/ic/progress-timeline.json`
- `docs/performance-gains/`

They record committed results in the group-addition unit, which stands as
recorded. The scoreboard row for the ECC2K-130 crossover is not a graph value.
Its number is withdrawn in `RESEARCH_ECC2K130_IC_FEASIBILITY.md`, with a
status note.

The study measures no speedup, so nothing goes in `benchmarks/`.

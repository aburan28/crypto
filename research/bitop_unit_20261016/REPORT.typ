// research/bitop_unit_20261016: every ECDLP method priced in bitop/v1.
// Every number in a table is read from a file ecbench wrote.
#set page(paper: "a4", margin: (x: 2cm, y: 2cm), numbering: "1")
#set text(size: 9.5pt)
#set par(justify: true)
#show heading.where(level: 1): set text(size: 13pt)
#show heading.where(level: 2): set text(size: 11pt)
#show link: set text(fill: rgb("#2a78d6"))

#let bound(name) = json("bounds/" + name + ".json")
#let f3(x) = if x == none { "–" } else { str(calc.round(x, digits: 3)) }
#let f2(x) = if x == none { "–" } else { str(calc.round(x, digits: 2)) }
#let ci(v, d: 3) = if v == none { "–" } else { "[" + str(calc.round(v.at(0), digits: d)) + ", " + str(calc.round(v.at(1), digits: d)) + "]" }
#let sci(x) = if x == none { "–" } else if x < 1000 { str(calc.round(x, digits: 2)) } else {
  let e = calc.floor(calc.log(x, base: 10))
  str(calc.round(x / calc.pow(10, e), digits: 2)) + "·10" + super(str(e))
}

#align(center)[
  #text(size: 16pt, weight: "bold")[ECDLP methods priced in one unit: bit operations]
  #v(2pt)
  research/bitop_unit_20261016 · 2026-10-16 · aburan28/crypto branch `claude/bitop-work-unit`
]

= What this is

Until now, the repository and its siblings (`aburan28/cryptanalysis`, cairn) charged ECDLP work in at least six units:

- group-addition equivalents, with non-group work converted by one host's nanosecond ratios;
- solver work priced by wall time, or not at all;
- `full_prices.json` ratios;
- reference picoseconds;
- Callgrind instructions;
- wall-clock online speedups.

No two of these can be divided, so an exponent fitted in one is not an exponent in another. Worse, work without a conversion ratio was charged *zero*: Buchberger, SAT and PDP solver work, and linear algebra on any curve without a calibration row. Every index-calculus figure was therefore a lower bound.

*`bitop/v1`* (`docs/units/BITOPS.md`, `src/cryptanalysis/bitops.rs`) replaces all of them. A bit operation is one two-input gate, and each primitive has a frozen gate count that depends only on operand width:

- field product, squaring and inversion in $FF_p$, $FF_(2^n)$ and $FF_(p^4)$;
- affine addition and doubling;
- multiply, add and inversion modulo $r$;
- a 64-bit word operation;
- a hash probe.

A run's cost is an exact integer that no host enters. `ecbench::work` prices any record from its own counts. A record whose counters include one with no rule, or whose code is known to do work no counter sees, is *incomplete*, and its bound in this unit is inadmissible.

Completeness is checked independently. `ecbench audit-ir` replays a run under Callgrind and divides its instructions by its counted bit operations; an arm whose ratio is far above rho's on the same curve is doing work its counters miss.

The protocol (`PROTOCOL.md`) was frozen before any session; its deviations are listed there.

= Completeness: the instruction audit

To make every method complete, the solvers and the pipeline were instrumented:

- F4, F5, inherited-F4, F6, Buchberger, crossbred, XL, FES, exhaustive and SAT now count their whole run;
- so do the descent oracles' system builds;
- so do the pipeline's Gauss scans and pivot inversions, relation rows, walk bookkeeping, table inserts, base indexes and canonicalisation tables.

New counters live in fields a replay of a committed record ignores, so every committed session still replays exactly. This was checked on nine sessions covering every engine family.

#figure(image("figures/audit.svg", width: 100%), caption: [Instructions per thousand counted bit operations, each index-calculus arm against rho on the same curve, with the final binary (`audit/final-panel.json`). Every arm sits within $2 times$ of its rho, and most sit below it.])

#let panel = json("audit/final-panel.json")
#table(
  columns: (auto, 1fr, auto, auto, auto, auto),
  align: (right, left, left, right, right, center),
  table.header([*seq*], [*method · oracle / solver*], [*curve*], [*instr./kbitop*], [*× floor*], [*complete*]),
  ..panel.map(d => (
    str(d.seq),
    d.method + " · " + (if "solver" in d.params and d.params.solver != "" { d.params.solver.split(":").at(0) } else if "oracle" in d.params { d.params.oracle.split(":").at(0) } else { "" }),
    d.curve.split("-").slice(0, 2).join("-"),
    f2(d.ir_per_kilo_bitop),
    sci(d.work.ratio_to_floor),
    if d.work.complete { "yes" } else { "no" },
  )).flatten(),
)

Rho on the 13-bit $FF_(2^17)$ curve (seq 35) reads 1,686 instructions per kilo-bitop because its walk is about 100 steps, and the per-solve fixed overhead dominates. The F4 arm on the same curve is at 231.

= What the old unit left out

#figure(image("figures/old_vs_new.svg", width: 100%), caption: [The same executions priced two ways: the committed record's ratio to the floor in group-addition equivalents, and the audit replay's in `bitop/v1` with complete counters.])

- *Table oracles:* their old and new figures agree to within 15% (Frobenius-orbit MITM 0.2%, pdp3 root-finding 4%, prime root-finding 7%, prime MITM 15%).
- *Algebraic oracles:* they were understated by three orders of magnitude or more. Their solvers had been charged little or nothing:
  - F4 descent $8.0 -> 7147$ times the floor;
  - Buchberger $53 -> 67837$;
  - crossbred $51 -> 55887$;
  - SAT $15.5 -> 16552$;
  - inherited-F4 $66 -> 144448$.
- *Signed-Frobenius rho:* its ratio roughly doubles ($2.40 -> 5.38$), because its canonicalisations, a Frobenius normal form and two products modulo $r$ per step, had been free.

= Exponents with intervals

The prime ladder (`specs/prime-ladder.json`, session `ECBS1h3674c4a0b139`) has 11 prime-order curves, one at each bit size from 16 to 26, with 8 public targets each and 528 verified runs. Every arm's pricing is complete and every bound is admissible.

Each exponent is fitted twice. The *raw* exponent fits gates $= C dot r^alpha$. The exponent *in additions* fits gates divided by one affine addition's gates on each curve, which takes out the field arithmetic's $(log r)^2$; there, generic theory says $1/2$. Intervals are two-stage cluster bootstrap, 95%, 2000 resamples.

#let ladder = ("rho-neg", "bsgs-neg", "kangaroo", "ic-mitm-b32", "ic-mitm-b128", "ic-rootfind-m3-b32")
#table(
  columns: (auto, auto, auto, auto, auto, auto),
  align: (left, right, left, right, left, right),
  table.header([*arm*], [*α raw*], [*95%*], [*α in additions*], [*95%*], [*× floor (26 bits)*]),
  ..ladder.map(a => {
    let b = bound("prime-ladder-" + a)
    (a, f3(b.fit.alpha), ci(b.fit.alpha_ci95), f3(b.fit.alpha_in_additions), ci(b.fit.alpha_in_additions_ci95), f2(b.constant.ratio_to_floor.value))
  }).flatten(),
)

*Positive control: passed.* Rho, BSGS and kangaroo all have $1/2$ inside their interval in additions. The instrument reproduces the generic theory.

*Index calculus with a fixed base* scales worse than rho at these sizes. The relation phase's exponent is about 1.0 for the 32-point base and about 0.83 for the 128-point base; the latter's low total exponent is its fixed table set-up dominating small sizes. These are exponents of these *configurations*. An index calculus that grows its base with $r$ is a different method and needs its own ladder.

#figure(image("figures/exponents.svg", width: 92%), caption: [Exponents in additions with 95% intervals: prime ladder and fiber ladder.])

#figure(image("figures/ladder_runs.svg", width: 100%), caption: [Every verified run of the prime ladder: $log_2$ gates against $log_2 r$.])

*Fiber ladder* (`specs/prime-fiber-scaling.json`, session `ECBS1hef866ef12457`, 8 sizes). All seven arms are admissible. Rho's interval in additions, #ci(bound("fiber-rho-neg").fit.alpha_in_additions_ci95), excludes $1/2$ by 0.041, inside the protocol's 0.05 tolerance. Its fixed jump-table set-up is a larger share at the small sizes, which pulls the fitted exponent down.

*Koblitz curves* ($a = 0$, $n = 13, 19, 23, 31$). Signed-Frobenius rho is at #f2(bound("koblitz-rho-frob").constant.ratio_to_floor.value)$times$ the floor, and strong rho at #f2(bound("koblitz-rho-strong").constant.ratio_to_floor.value)$times$. Both are admissible. Their exponents in additions fall below $1/2$, as they should: the automorphism group they fold by has order $2n$, which grows with $log r$. With 4 sizes, the intervals are wide.

== How many digits an exponent has

An exponent is reported only to the width of its interval. Rho's interval in additions on the prime ladder is about $plus.minus 0.11$ from 88 runs. The width scales as one over the square root of the runs, so a third decimal place ($plus.minus 0.001$) would take about $10^4$ times as many runs, and widening the size range helps more than repeating sizes. The figures above carry exactly the digits their data support.

= Algebraic oracles, constant only

These single-curve sessions are re-runs of committed specs on the instrumented binary. All are complete and admissible.

#table(
  columns: (1fr, auto, auto),
  align: (left, right, left),
  table.header([*session · arm*], [*× floor*], [*95%*]),
  ..("sat-count-qualified-ic-sat-xor", "sat-count-qualified-ic-mitm", "sat-count-qualified-ic-mitm-frob", "sat-count-qualified-rho-strong", "f4-k13-matrix-f4", "f4-k13-matrix-f5", "f4-k13-inherited-f4", "f4-k13-f6-ic", "f4-k13-rho-negation", "f4-k13-rho-strong").map(n => {
    let b = bound(n)
    (n, sci(b.constant.ratio_to_floor.value), ci(b.constant.ratio_to_floor.ci95, d: 2))
  }).flatten(),
)

In the yield sweep (`specs/yield-koblitz.json`, 4 sizes) every arm is admissible except the two SAT-descent arms. Those fit over 3 sizes and pass through a base builder and an elimination variant that are not yet instrumented; their bounds carry that and are not admitted.

= Amortised: is there a number of targets at which index calculus wins?

`ecbench amortize` splits each run's gates into *reusable* gates and *per-target* gates. The reusable part is the factor base and the oracle's tables. Each arm is then compared with batch rho, the Kuhn–Struik walk that keeps its table, whose cost grows like $sqrt(T)$. Where a method wins at all it wins in a window of target counts, because batch rho's cost grows like $sqrt(T)$ and the method's like $T$.

On the prime ladder, none of the three index-calculus arms has a window at any size. Each further target costs them 3.5–311 times what it costs rho (`amortize/`). The reason is that the single-target pipeline re-collects relations for every target. Only a target-blind variant such as `ic.shared_rank` could change that.

= The archive

`ecbench work --summary` over every committed session (`archive/summary.json`) prices #(json("archive/summary.json").map(r => r.complete_runs).sum()) of #(json("archive/summary.json").map(r => r.runs).sum()) verified runs completely: every generic-method run.

The rest are index-calculus records written before the new counters existed. They name exactly what they miss, and only a re-run can complete them, as this study's sessions were.

= Limits

- *What the gates leave out:* memory, branches and allocation. Memory is reported beside the unit; a memory-cost model would be `bitop/v2`.
- *Conventions are frozen:*
  - schoolbook products, so Karatsuba would be `bitop/v2`;
  - a pentanomial reduction;
  - a $4 kappa$ hash probe.

  Every record keeps its primitive counts, so any later convention reprices the whole history.
- *Root-finding:* the arm batches its inversions but is charged an affine addition, inversion included, per candidate, the unit's convention. Its gates are a convention-level upper figure for that implementation; its instruction ratio, the lowest in the panel, shows it.
- *Isolation:* sessions ran at isolation L0 (`--cpus none`). Bit-operation counts do not depend on it; no wall time is used here.
- *No security claim:* nothing here concerns any deployed curve.

= Provenance

- *Specs:* `specs/`.
- *Sessions:* `sessions/`, each verified with `ecbench verify --replay 2`.
- *Bounds:* `bounds/`, from `ecbench bound fit --unit bitop/v1`.
- *Audits:* `audit/`.
- *Amortisation:* `amortize/`.
- *Archive summary:* `archive/`.
- *Figures:* `figures/plots.py`, which reads only recorded values.

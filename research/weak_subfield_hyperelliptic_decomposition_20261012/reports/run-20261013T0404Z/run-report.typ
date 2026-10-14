#set page(paper: "a4", margin: 1.8cm, numbering: "1")
#set text(size: 9.5pt)
#show table: set text(size: 8pt)
= Run report: weak-curve laboratory, run 20261013T0404Z

Descriptive record of one run. Interpretation is in the study `README.md` / `REPORT.pdf`.

== What ran
- Command: `rustc --edition 2021 -O weak_curve_lab.rs && ./weak_curve_lab all && ./weak_curve_lab e4x` (experiments E5, E1, E2, E3, E3b, E4 from `all`, then E4x from the second command; `all` no longer includes E4x, PROTOCOL.md amendment 1).
- Corrected re-run of `20261013T0055Z`, which recorded E4x twice (`all` included it). Records 1-116 and the three E4x records are identical to that run's (checked line by line). Spec: `PROTOCOL.md`, frozen before the run, plus dated amendments. E4x (sextic with `a, b >= 1`) was added after E4 was first seen and is supplementary.
- Environment, base commit, source SHA-256, rustc version, host, load and wall time: `results/run-20261013T0404Z.env.txt` (reproduced below).
- Isolation: none. Shared host under heavy load; no timing is reported as a result. All recorded values are exact integer counts and are deterministic (fixed seed 20261012).
- Raw records: `results/run-20261013T0404Z.jsonl`, one JSON object per line, field `exp` names the experiment.

#raw(read("../../results/run-20261013T0404Z.env.txt"), block: true)

== Runs requested and completed
#table(columns: 4, stroke: 0.4pt,
  [*Experiment*], [*Requested*], [*Completed*], [*Notes*],
  [E5 paper data], [2 examples], [2], [],
  [E1 collapse + control], [5 fields x (3 weak + 3 random)], [30], [no subgroup hit the cap; samplers draw nonzero coefficients (amendment 2)],
  [E2 lines], [4 fields x 13 lines], [52], [],
  [E3 census], [6 fields x 4 curves], [23], [`x^5+c` not run at (5,4): singular in characteristic 5],
  [E3b multi-form, post-hoc], [6 curve/field pairs], [6], [amendment 3; union capped at 1.2e7 elements for (11,5) and (13,4); reported as lower bounds],
  [E4 Kani], [3 fields], [3], [],
  [E4x Kani, supplementary], [3 fields], [3], [not in the frozen protocol; run once by `e4x`],
)

== Recorded values
All values are *measured* exact counts unless marked derived. Unit: number of group elements (`H`, `J`), number of pairs (`candidates`), number of lines.

#figure(image("e1_span_fraction.svg", width: 100%), caption: [E1, run 20261013T0404Z. Each marker is one curve (n = 3 per field and kind). Vertical axis: log |H| / log \#J(K), derived as a ratio of two measured counts. Dashed: the value gcd(2,n)/n. No uncertainty: counts are exact.])

#figure(image("e2_lines.svg", width: 100%), caption: [E2, run 20261013T0404Z. One weak curve per field, 13 factor-base lines each. Top: measured number of relation candidates (symlog axis; Moebius lines use the sextic model and are not like-for-like with the others). Bottom: log |H| / log \#J(K).])

#figure(image("../mechanism.svg", width: 100%), caption: [Diagram of what was run and where the measured subgroups sit.])

== Controls and cross-checks recorded in the run
- `orders_agree`: \#J(K) from point counts over K and K^2 against det(1 - pi^n) from the F_q L-polynomial, on every weak curve and every Kani curve.
- `selftest`: \#J(K) D = 0 and associativity on 20 random triples per E1 curve.
- `H_divides_sub`, `coords_in_subfield` (weak curves); `H_divides_J` (random curves).
- `hom_ok` / `hom_tested` (E4): homomorphism property of the pullback on random pairs that avoid the fibre `u = r^2`, where the map is not implemented.
The per-record values are in the JSONL; the study report tabulates them.

== Missing data
- Union subgroup order for `x^5+c` at (11,5) and `x^5+cx` at (13,4): capped, lower bound recorded.
- Span of the second weak line found by the E3 census at (7,4): not computed.

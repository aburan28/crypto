# Different-j full-image factor-base gate: result and decision

**Status:** the preregistered one-target full-image construction,
three-round paired relation solve, exact support census and independent
replay are complete on the frozen `F_(67²)`, order-61 pair. The
[run report](reports/run-20261015T1916Z/REPORT.md),
[visual PDF](reports/run-20261015T1916Z/REPORT.pdf),
[V1 protocol](PROTOCOL_V1.md), [hash-pinned V2 spec](RUN_SPEC_V2.json),
[raw records](RUN_V2.jsonl), [separate replay](REPLAY_V2.jsonl),
[native phase summary](SUMMARY_V2.json) and
[exact subgroup support](SUPPORT_V2.json) carry the evidence.
The run is stacked on PR #1920 and is not an ecbench frontier claim.

| Requested gate | Status | Evidence and limit |
| --- | --- | --- |
| Full rational image and deterministic base | Verified for the frozen pair | 4,575 source points, 1,525 image points, three preimages each; eight projected columns selected from 21 sorted candidates; eight source witnesses replayed. |
| Paired relation collection and scalar recovery | Verified for this public target | Four relation arms reached rank nine in 3/3 rounds; all 18 relation and rho arm-round rows recovered scalar 41 and passed separate scalar replay. |
| Exact support attribution | Verified for the four fixed bases | Each base covers 56 of 60 nonzero sums; full-image has 52 common sums and four unique sums relative to both projected image and destination direct. Source and projected image support sets coincide. |
| Complete local timer-interval cost | Measured at shared-host L0 | Full-image median current-thread CPU 3,962.083 µs versus direct rho 58.084 µs and destination direct 134.584 µs; paired full/rho median ratio 70.281 and full/direct 29.959. Full-image setup median share 98.532%. |
| Portable method price and scaling | Open | Common input/model/startup work is outside the custom timer; field/map/lookup work is not converted to ecbench's counted `S`; isolation, independent-host cost replay and increasing-order tests are absent. |

The structural observation is exact for this factor-base selection:
equal 56-point support counts did **not** mean equal support membership.
Full-image versus projected-image had four sums unique to each, and
full-image versus destination-direct also had four unique to each.
The [Lean finite-support proof map](PROOF_MAP.md) checks the exact
residue sets and 52/4/4 differences for the frozen scalar-coordinate
lists with no axioms. The bridge from the listed curve points to those
scalar coordinates and the concrete Vélu map remains a separate
curve-level formalization obligation; native arithmetic replay checks
the bounded instance.
That explains why relation trajectories may differ even when support
cardinality does not. The measured full-image relation `(trials, hits)`
were `(11,10)`, `(23,19)`, `(11,11)`; the destination-direct arm had
`(15,14)`, `(11,11)`, `(21,19)`. The three rounds do not isolate a
stable relation-yield change.

The measured cold total was higher because the selected full image was
constructed on every target run. Its setup used 9,343 group additions,
390,029 extension-field multiplications and 13,819 inversions per
round, compared with 136 additions, 7,434 multiplications and 299
inversions for the direct destination base. The fully charged local
timer ratio to same-target rho was 73.826, 70.281 and 22.743 by round;
the median is **70.281**, calculated from paired round ratios rather
than from the displayed arm medians. This is a bounded engineering and
coordinate-support result; it supplies no exponent improvement.

![Full-image cold CPU and exact support comparison](reports/run-20261015T1916Z/GRAPH.svg)

The next discriminating test is a **separately frozen 60-target reuse
panel**: charge the 1,525-point image construction once, compare cold
and truly amortized totals with direct source BSGS and rho, and replay
every scalar. Its [protocol](NEXT_GATE.md) keeps the still-open counted
`S` and independent-host obligations explicit. The canonical frontier
and progress charts in `docs/bounds/frontier.json`,
`docs/index-calculus-scoreboard.html` and
`docs/ic/progress-timeline.json` were checked and left unchanged because
this run adds no admitted counted-`S` crossover. The study graph above
is the new quantitative record.

The initial `serde_json::json!` macro recursion-limit build failure,
corrected source/binary hashes, zero-test example harness and exact
commands are retained in [BUILD_V2.md](BUILD_V2.md). The
[validation receipt](VALIDATION.md) records the controlled release
library suite: 4,083 passed, zero failed and 117 ignored in 957.05 s;
the three touched example harnesses compiled, and the finite Lean
support module built with no axioms. The run script preserved 18
individual stdout, stderr and exit-code triples in `raw-v2/`; all
18 exited zero.

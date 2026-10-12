#set page(paper: "a4", margin: 1.8cm)
#set text(font: "DejaVu Sans", size: 9pt)
#set par(justify: false)
#set table(stroke: 0.4pt + rgb("#c3c2b7"), inset: 4pt)
#show table: set text(7.6pt)
#text(13pt, weight: "bold")[Run report: run-20261011T2330Z-gls-validate]

*What ran.* `./target/release/examples/gls_loop_validate --pbits 127 --points 24 --seed 20261012 --discs -619,-23,-71,-431,-3019 --json research/isogeny_loop_glv_20261012/gls_validate.json` (in `isogeny_algos/`). Source: working tree that became checkpoint commit 54a16f6d4 on branch `claude/isogeny-loop-glv-sweeper`. rustc 1.97.0, release profile, Linux 6.18 x86_64 cloud container, not isolated. Start 2026-10-11T23:30Z (`../../validate_run_start_utc.txt`); wall time 2.8 s (`../../gls_validate.log`).

*Replay.* On 2026-10-12 the same command was re-run on the final source, after helper functions moved into `src/gls_loop.rs` and the clippy-only edits. The Markdown table was byte-identical, and the JSON differed only in its recorded `command` (output path).

*Arms.* One arm per discriminant, D in {−619, −23, −71, −431, −3019}. Each runs these steps:
- certified cheapest loop;
- instance search (4p = t² + |D|v², r = \#E′(F#sub[p²])/h prime, p ≡ 3 mod 4);
- H#sub[D] by the analytic method, roots mod p, curve, and a twist if needed;
- crater walks;
- an F#sub[p] check on E(F#sub[p])[r#sub[E]];
- the ψ and twisted-loop checks on 24 points of order r on E′.

*Runs requested / completed.* 5 / 5, exit 0, verdict ALL PASS. No failures.

*Recorded values* (`../../gls_validate.md`, `../../gls_validate.json`):

#table(columns: 7,
  [*D*], [*steps*], [*H#sub[D] deg / bits / worst err*], [*log#sub[2] r*], [*ψ*], [*loops act as*], [*verdict*],
  [−619], [5,5,7], [5 / 303 / 2#super[−90]], [251], [+μ 24/24], [−λ(s); −λ(−s)], [PASS],
  [−23], [2,3], [3 / 150 / 2#super[−86]], [251], [+μ 24/24], [−λ(−s); −λ(s)], [PASS],
  [−71], [2,3,3], [7 / 238 / 2#super[−92]], [248], [+μ 24/24], [−λ(−s); −λ(s)], [PASS],
  [−431], [2,2,3,3,3], [21 / 631 / 2#super[−118]], [252], [+μ 24/24], [+λ(−s); +λ(s)], [PASS],
  [−3019], [5,5,5,13], [7 / 523 / 2#super[−88]], [253], [+μ 24/24], [−λ(−s); +λ(s)], [PASS],
)

#figure(image("../../figures/status.svg", width: 92%), caption: [What this run validates (green, top right) within the round-3 pipeline. Diagram, not data: each row has a single pass/fail value, so there are no comparable quantities to plot.])

*Observations only.* Interpretation: `../../README.md` §2.

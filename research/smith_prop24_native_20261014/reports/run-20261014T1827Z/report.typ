// Run report (descriptive). Every number is read from the recorded files.
#set page(paper: "a4", margin: 1.6cm, numbering: "1")
#set text(size: 9pt)
#set par(justify: true)
#let d = "/data/"
#let man = json(d + "run_manifest.json")
#let q1 = json(d + "random_p2147483647_summary.json")
#let c7 = json(d + "census_p7_summary.json")
#let c11 = json(d + "census_p11_summary.json")
#let r7 = json(d + "census_p7_raw.json")
#let r11 = json(d + "census_p11_raw.json")
#let ctl = json(d + "control_example1.json")
#let f(x, n: 4) = str(calc.round(x, digits: n))
#let tbl(..args) = text(size: 7.5pt, table(inset: 3pt, stroke: 0.4pt + luma(170), ..args))

#align(center)[#text(size: 14pt, weight: "bold")[Smith's Propositions 2 and 4, natively: run report] \
Run `run-20261014T1827Z` (UTC start) · study `research/smith_prop24_native_20261014` · descriptive, no verdicts]

*What ran.* `examples/smith_prop24.rs` built with `cargo build --release --example smith_prop24` at commit
#raw(man.at(0).git_head.slice(0, 9)) (the protocol freeze commit; `git status` of `src/` and `examples/` recorded as
#if man.all(r => r.git_dirty_paths == "") [clean] else [NOT clean] at the end of all #man.len() runs). `PROTOCOL.md` SHA-256 #raw(man.at(0).protocol_sha256.slice(0, 16))…, the same in every
manifest entry. Host: 4-vCPU Intel Xeon \@ 2.80 GHz cloud container, 15 GB RAM, Linux 6.18, #raw(man.at(0).rustc), rayon with
#man.at(0).threads threads, shared with another agent's builds, no isolation. Wall times are incidental; the results are exact
computations and counts, not timings.

#tbl(columns: 4,
  [*run*], [*command (after the binary)*], [*seconds*], [*status*],
  ..man.map(r => ([#r.run], [#raw(r.command.slice(1).join(" "))], [#f(r.seconds, n: 1)], [completed])).flatten())

*Arms and seeds.* Control: Smith's Example 1 over $FF_37$ and the same curve in a degree-8 model. Question 1:
$N = #q1.n_curves$ random squarefree degree-8 models over $p = 2^31 - 1$, seed #q1.seed (`SplitMix64`, curve $i$ from
`seed XOR i·0xD1B54A32D192ED03`). Question 2: every monic octic with zero $x^7$ coefficient over $FF_7$ (#r7.octics_enumerated codes)
and $FF_11$ (#r11.octics_enumerated codes), each squarefree one as a curve and its twist by #r7.twist_constant (p = 7) and
#r11.twist_constant (p = 11); 199 permutations of the independence null with seed #c7.null.seed. Requested runs: 4. Completed: #man.len().
Failed or timed out: none. After the runs, `smith_prop24 resummarize` rebuilt every summary byte for byte from the rows,
and `verdict` and `describe` (post hoc, protocol amendment 1) read only the recorded files.

*Measured values* (unit in each row; intervals are those the summary recorded).

#tbl(columns: 6,
  [*Q1 fraction*], [*k / n*], [*value*], [*98.75% Clopper–Pearson*], [*prediction*], [*Smith Ex. 2*],
  ..q1.fractions.map(x => ([#x.name], [#x.k / #x.n], [#f(x.value, n: 5)], [[#f(x.cp_lo, n: 5), #f(x.cp_hi, n: 5)]],
    [#f(x.prediction, n: 5)], [#f(x.smith_example2, n: 5)])).flatten())

#tbl(columns: 9,
  [*p*], [*curves*], [*classes*], [*direct full\_success D*], [*class-level C*], [*ceiling U*], [*r*],
  [*null C (2.5–97.5%)*], [*tractable / full classes*],
  ..(c7, c11).map(s => ([#s.p], [#s.curves], [#s.classes], [#f(s.direct_full)], [#f(s.class_full)], [#f(s.class_tractable)],
    [#f(s.r, n: 3)], [#f(s.null.class_fraction_lo)–#f(s.null.class_fraction_hi)],
    [#s.classes_with_tractable / #s.classes_with_full])).flatten())

*Instrument diagnostics* (counts over all subgroups or maps of the run; the protocol requires the fail, mismatch,
disagreement and error columns to be 0).

#let dg = (("Q1", q1.diag), ("p = 7", c7.diag), ("p = 11", c11.diag))
#tbl(columns: 12,
  [*run*], [*maps*], [*rank M < 4*], [*disc = 0*], [*not normal form*], [*deg g < 3*], [*construction failed*],
  [*pair fail*], [*eq. (7) vs norm*], [*Lemma 4 fail*], [*λ disagree*], [*X agree / disagree / n.a. / error*],
  ..dg.map(((n, x)) => ([#n], [#x.maps], [#x.m_rank_deficient], [#x.disc_zero], [#x.non_normal_form], [#x.degree_drop],
    [#x.construction_failed], [#x.pair_condition_fail], [#x.eval_mismatch], [#x.lemma4_fail], [#x.lambda_disagree],
    [#x.x_agree / #x.x_disagree / #x.x_not_applicable / #x.x_count_error])).flatten())

#let yn(b) = if b [yes] else [NO]
*Controls.* Example 1 passed: #yn(ctl.passed) (one rational tractable subgroup, Smith's $N, D$ among the computed maps,
$delta_0, delta_2, delta_4$, $64 s = delta_1^2$, $f_0, f_1, f_2$ and $\#H(FF_37) = #ctl.h_points$ all as printed by Smith, Prop. 4
rational in both models). Frozen 2026-10-13 counts reproduced: p = 7 #yn(c7.control_frozen_counts_match), p = 11
#yn(c11.control_frozen_counts_match). Class closure of "has a rational tractable subgroup" equals its direct fraction:
p = 7 #yn(c7.control_tractable_closure_equal), p = 11 #yn(c11.control_tractable_closure_equal).

#figure(image("/figures/smith_prop24_results.svg", width: 100%),
  caption: [Recorded values of this run. (a) Question 1, observed minus predicted, 98.75% Clopper–Pearson intervals.
  (b) Question 2, exhaustive: D, C, U per prime. Source: `figures/plot.py` from `data/*.json`.])

*Raw records.* `data/random_p2147483647_rows.csv.gz` (one row per curve), `data/census_p{7,11}_classes.csv` (one row
per class), `data/census_p{7,11}_types.csv`, `data/census_p{7,11}_raw.json`, the summaries
`data/*_summary.json`, `data/control_example1.json`, `data/run_manifest.json`. Interpretation and verdicts are in
`REPORT.pdf`, not here.

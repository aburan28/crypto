#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9.2pt)
#set table(stroke: 0.4pt)
#show table: set text(size: 8pt)
#show figure.caption: it => it.body

#let S = json("summary_edwards.json")
#let sci(x, d: 2) = {
  if x == none { return "—" }
  let v = float(x)
  if v == 0 { return "0" }
  let e = calc.floor(calc.log(calc.abs(v), base: 10))
  let m = v / calc.pow(10.0, e)
  [#calc.round(m, digits: d)·10#super[#e]]
}
#let r2(x, d: 2) = if x == none { "—" } else { str(calc.round(float(x), digits: d)) }

= $k = 5$, 2-torsion Edwards route, end to end with F4 trace replay

#box(fill: rgb("#fff4e0"), inset: 6pt, radius: 3pt)[*Provisional.* Built from the records present when it was rendered: DLP arms without a final row are shown as incomplete. The protocol's remaining arms are still running.]

2026-10-15 · `research/joux_vitse_n1_20261015/` · protocol `PROTOCOL_EDWARDS.md` (frozen at `3a8599b3f`), deviations `data/edwards/DEVIATION.md` · code `src/cryptanalysis/joux_vitse_edwards.rs`, `examples/joux_vitse_edwards.rs`, summary `examples/joux_vitse_edwards_summary.rs` · construction unchanged from `jv_quintic_edwards.rs` (RESEARCH_K5_TORSION_JOUX_VITSE.md). Checks failed so far: *#S.checks_failed*.

== The test: trace replay against the committed solve, same residuals

#table(columns: 9, align: (left,) + (right,) * 8,
  [run], [$|F|$], [committed solve (mean)], [replay (mean)], [factor mean], [min], [max], [replay div.], [checks failed],
  ..S.cost.map(c => (
    c.file, str(c.base_points), sci(c.full_mean), sci(c.replay_mean), r2(c.factor_mean), r2(c.factor_min), r2(c.factor_max), str(c.replay_mismatches), str(c.checks_failed),
  )).flatten()
)

On the committed curves ($p = 271, 521$, the bases of RESEARCH_K5 §3), the committed test now costs $5.85 dot 10^7$ (the note's $7.3 dot 10^7$ was measured on the 2026-09-29 engine) and the replay $3.3$–$3.4 dot 10^7$: a factor *1.74–1.77*, below the $3$–$10 times$ the note priced. Every replayed residual's decompositions equal the committed solve's and the meet-in-the-middle oracle's.

#figure(image("figures/edwards_replay.svg", width: 86%), caption: [Cost of one test per residual: committed solve and replay (mean, min–max), with the note's figure.])

== Complete DLPs, paired with rho on the same target

#table(columns: 10, align: (left,) + (right,) * 9,
  [run], [$log_2 n$], [$|F|$], [trials], [rel.], [trials/rel.], [÷ census], [IC online], [IC / rho], [IC / exp. walk],
  ..S.dlp.map(d => if d.complete {(
    d.file, r2(d.log2_n, d: 1), str(int(d.base_points)), str(int(d.trials)), str(d.relations), r2(d.trials_per_relation, d: 0),
    r2(d.trials_per_relation_over_census) + " [" + r2(d.trials_per_relation_over_census_95.at(0)) + ", " + r2(d.trials_per_relation_over_census_95.at(1)) + "]",
    sci(d.online_muls), sci(d.ratio_online_ic_over_rho), sci(d.ratio_online_ic_over_expected_rho_walk),
  )} else {(
    d.file, "", str(int(d.base_points)), "running", str(int(d.relations_recorded)) + " so far", "—", "census " + r2(d.census_trials_per_relation, d: 0), "—", "—", "—",
  )}).flatten()
)

The census column is the measured trials per relation over the census expectation $n \/ (8 binom(|F|, 4))$ (RESEARCH_K5 §3), with an exact Poisson 95 % interval on the relation count; it is an expectation, not a bound. At $p = 31$ the complete DLP recovered the planted logarithm (verified $k G = Q$; rho on the same target verified too) at *#sci(S.dlp.at(0).ratio_online_ic_over_rho) × rho*; the plain Weierstrass route at the same $p$ measured $6.8 dot 10^7$. About 3 % of residuals ($approx 1\/p$) give a consistent system; each pays the replay and the full solve, so this run's mean trial is #sci(S.dlp.at(0).c_trial) against the replay's $3.4 dot 10^7$. `dlp_p41_s1` drew a six-point base (census $approx 2.4 dot 10^5$ trials per relation) and is expected to reach its 10-hour cap.

== Where it would meet rho (extrapolation)

On the measured replay cost (#sci(S.model_extrapolation.c_trial_replay)), the census rate $1\/(192 p)$, $|F| + 1 approx p\/4$ relations and rho with the negation map at $c_"add" = 452$: parity at $n approx 2^(#r2(S.model_extrapolation.log2_n_star, d: 1))$ ($p approx 2^(#r2(S.model_extrapolation.log2_p_star, d: 1))$), against the note's $2^233$ without replay. At the measured $p = 31$ ratio to the census (#r2(S.model_extrapolation.measured_ratio_to_census), interval containing 1) it would read $2^(#r2(S.model_extrapolation.log2_n_star_at_measured_ratio, d: 1))$; that is a sensitivity, not the estimate. Either way the route stays $10^5$–$10^6$ times rho at every size this repository can run.

#figure(image("figures/edwards_ic_over_rho.svg", width: 86%), caption: [IC / rho: complete DLPs (Edwards, circles; plain Weierstrass, squares) and the two models (dashed; extrapolation).])

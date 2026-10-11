// lit2026 W3 results report. Every number is read from the run's JSON
// (written by examples/lit2026_w3_ggm_audit.rs); nothing is computed here
// beyond formatting.
#let run = sys.inputs.at("run", default: "reports/run")
#let t = json(run + "/table.json")
#let fit = json(run + "/fit.json")
#let meas = json(run + "/measure.json")
#let meta = json(run + "/meta.json")
#set page(paper: "a4", margin: 1.6cm, numbering: "1")
#set text(size: 9pt)
#set par(justify: true)
#show heading: set block(above: 1.1em, below: 0.6em)

#let g(x, d: 3) = if x == none { "—" } else if type(x) == str { x } else {
  let a = calc.abs(x)
  if a == 0 { "0" } else if a >= 1e6 or a < 1e-3 {
    let e = calc.floor(calc.log(a, base: 10))
    str(calc.round(x / calc.pow(10.0, e), digits: 2)) + "e" + str(e)
  } else if a >= 100 { str(calc.round(x)) } else { str(calc.round(x, digits: d)) }
}

= lit2026 W3: IC and precompute configurations against generic lower bounds

*Run* #raw(meta.run_id) · *commit* #raw(meta.commit) · *host* #meta.host · *date* #meta.date.
Papers: Corrigan-Gibbs–Henzinger–Wu, _The Structured Generic-Group Model_, ePrint 2026/384 (EUROCRYPT 2026);
Dinur–Keller–Marmor, arXiv 2505.00894 (STOC 2026); Akshima–Kiltz–Nishaburi–Nooripoor–Wiedijk, ePrint 2026/2253 (TCC 2026).
Producer: `examples/lit2026_w3_ggm_audit.rs` (`table`, `fit`, `measure`). Plot: `plot.py` (layout only).
Scope: public synthetic toy curves; no key-recovery, asymptotic or deployed-curve claim; multi-target results are secondary diagnostics.

== Bounds used

- *SGGM, Theorem 3.2.* Prime order $r$; $delta$ = fraction of labels constrained by $star$ (Def. 3.1: inputs or output of a defined $star$). For every $T$-query algorithm, $"Adv" <= delta(3T+2) + (3T+1)^2 / r + 1/r$, so constant success needs $T = Omega(min(sqrt(r), 1/delta))$ group-operation queries; $star$ evaluations are free. Appendix B (Thm B.1): with $S$ advice bits, $"Adv" <= tilde(O)(S T^2 / r + delta T)$.
- *Mapping to IC.* Labels = points of the order-$r$ subgroup. An $m$-term decomposition oracle over factor base $F$ is a partial $star$ defined on pairs of $j$-sums ($j <= m$) whose sum is an $(i+j)$-sum; the constrained set is $F$ and every $j$-sum of $F$, so $delta <= c sum_(j=1)^m binom(|F|+j-1, j) \/ r$ ($c$ = torsion classes a relation may absorb; we take $c = h$). We use $delta = max("that count", "measured relations/trials")$, the weakest bound any record admits. The paper's copy we read (UT Austin PDF of the EUROCRYPT version, 2026-02-24) has no elliptic-curve section; the abstract's EC claim could not be checked because eprint.iacr.org sits behind a bot challenge.
- *Automorphisms.* Thm 3.2 has none. We report $min(sqrt(r\/A), 1\/delta)$ as asked and $min(sqrt(r\/A), 1\/(A delta))$ (each query read as $A$ labels; our extension, the weaker).
- *DKM Theorem 1.1.* Non-adaptive $(S,T)$ (pre-translated queries fixed by the advice, not by the challenge): $"Adv" <= 3T^2\/N + 4 ln 2 dot S T\/N$, so $S T gt.tilde N$ below $T = sqrt(N)$. Adaptive algorithms reach $S T^2 approx N$ (Corrigan-Gibbs–Kogan; Bernstein–Lange).
- *Multi-instance.* $k$ generic logarithms cost $Theta(sqrt(k r))$ in total (Yun 2015; 2026/2253 abstract: $1/2(log p + log k)$ bits for multi-instance CDH/SDH), so $d ln("total") \/ d ln k -> 1\/2$; Kuhn–Struik's law gives the finite-$k$ slope.

== Result 1: no measured cost is below its bound

#let fams = t.rows.map(r => r.family).dedup()
#let fam-rows = ("prime", "char2", "koblitz_round5", "koblitz", "koblitz_batch", "koblitz_compact_orbit")
#table(
  columns: 6, align: (left, right, right, right, right, right), stroke: 0.4pt + luma(180),
  table.header([family], [rows], [$delta_"upper" < 1$], [min measured / bound], [min measured / bound, $delta_"upper"<1$], [rows below a bound]),
  ..fam-rows.map(f => {
    let rs = t.rows.filter(r => r.family == f)
    let nv = rs.filter(r => r.delta_counting_upper < 1)
    let below = rs.filter(r => r.flags.any(x => x.starts-with("below"))).len()
    (f, str(rs.len()), str(nv.len()), g(calc.min(..rs.map(r => r.ratio_measured_over_theta))),
     if nv.len() > 0 { g(calc.min(..nv.map(r => r.ratio_measured_over_theta))) } else { "—" }, str(below))
  }).flatten()
)
Units: GAE = group-addition equivalents, whole process ($S dot sqrt(r) dot k$, $S$ from `docs/ic/leaderboard.json` and `docs/ic/boundary_targets.json`); compact-orbit rows are probes (S3 solve + root-table lookup), not group operations. Of #t.rows.len() rows, #t.rows.filter(r => r.delta_counting_upper >= 1).len() have $delta_"upper" = 1$: their factor base's $m$-sums cover the subgroup, so the bound is $min(sqrt(r\/A), 1) = 1$ and says nothing.

#figure(image(run + "/graph.svg", width: 100%), caption: [Recorded values from run #raw(meta.run_id). (a) every leaderboard/round-5 configuration, bound vs measured GAE; hollow = vacuous ($delta_"upper"=1$). (b) compact-orbit ladder ratios. (c) per-target multi-target rho cost relative to $k=1$, committed records (16/16/8 batches) and the fresh run (#meas.batches_per_k batches per $k$). (d) multiplier closure: canonical entries × mean first-hit probes over $r\/A$.])

=== Rows where the bound is not vacuous ($delta_"upper" < 1$)

#let nvrows = t.rows.filter(r => r.delta_counting_upper < 1)
#text(size: 7pt, table(
  columns: 10, stroke: 0.3pt + luma(190), align: (left,) + (right,) * 9,
  table.header([configuration], [log2 r], [A], [m], [|F|], [$delta$ meas.], [$delta_"upper"$], [bound], [measured GAE], [ratio]),
  ..nvrows.map(r => (r.id, str(calc.round(r.log2_r, digits: 2)), str(r.aut), str(r.m), str(r.factor_base_points),
    g(r.delta_measured), g(r.delta_counting_upper), g(r.bound_theta), g(r.measured), g(r.ratio_measured_over_theta))).flatten()
))

=== Measured relation rate above the counting ceiling

Relations per trial is a mean number of representations; for uniform targets its expectation is at most the count. #t.rows.filter(r => r.rate_above_ceiling_test != none).len() rows exceed it; tested as a Poisson excess over all trials (one-sided, Bonferroni level 0.05/100):
#text(size: 7.5pt, table(
  columns: 5, stroke: 0.3pt + luma(190), align: (left, right, right, right, right),
  table.header([row], [observed relations], [expected at ceiling], [z], [p]),
  ..t.rows.filter(r => r.rate_above_ceiling_test != none).map(r => (r.id,
    g(r.rate_above_ceiling_test.observed_relations), g(r.rate_above_ceiling_test.expected_at_ceiling),
    g(r.rate_above_ceiling_test.z), g(r.rate_above_ceiling_test.p_one_sided))).flatten()
))
None is significant after correction; the smallest $p$ is the shared `fp20` random-target trial stream (four variants report the same trials).

== Result 2: the compact-orbit ladder (n37–n83)

#let cr = t.rows.filter(r => r.family == "koblitz_compact_orbit" and r.id.ends-with("target_probes_online"))
#text(size: 7.5pt, table(
  columns: 10, stroke: 0.3pt + luma(190), align: (left,) + (right,) * 9,
  table.header([rung], [log2 r], [B], [K], [$delta_"upper"$], [target probes], [over $sqrt(r\/A)$], [log2(S·T²) − log2(r/A)], [S·T/(r/A) frozen / rank mean], [frozen / rank mean (P)]),
  ..cr.map(r => (r.id.split("/").at(0), str(calc.round(r.log2_r, digits: 2)), str(r.factor_base_points), str(r.orbit_columns),
    g(r.delta_counting_upper), g(r.measured), g(r.probes_over_sqrt_r_over_aut), str(calc.round(r.log2_s_t2_minus_log2_r_over_a, digits: 1)),
    g(r.s_t_over_r_over_a) + " / " + g(r.s_t_over_r_over_a_at_rank_mean),
    g(r.target_over_rank_mean) + " (" + g(r.p_exponential_at_most_observed) + ")")).flatten()
))
$S = K^2 n$ regular states of the S3 root index (advice); $T$ = probes. $delta_"upper" = 1$ on every rung ($B^4\/24 > r$), so the SGGM bound is vacuous: the cost lies in finding a decomposition, which the model treats as a free $star$ evaluation. Against preprocessing lines the online phase is consistent: $log_2(S T^2) - log_2(r\/A)$ is #g(calc.min(..cr.map(r => r.log2_s_t2_minus_log2_r_over_a)), d: 1) to #g(calc.max(..cr.map(r => r.log2_s_t2_minus_log2_r_over_a)), d: 1) with the frozen target's probes, and at the rank-stage mean probe count $S T \/ (r\/A)$ is #g(calc.min(..cr.map(r => r.s_t_over_r_over_a_at_rank_mean))) to #g(calc.max(..cr.map(r => r.s_t_over_r_over_a_at_rank_mean))), i.e. the target scan sits on the *non-adaptive* line $S T approx N$.

== Result 3: precompute arms and adaptivity

#let mc = t.precompute_multiplier_closure
Multiplier closure (#mc.len() arms, cells C23/C37/C37b): canonical entries × mean first-hit probes / $(r\/A)$ ranges #g(calc.min(..mc.map(p => p.canonical_times_t_over_r_over_a))) – #g(calc.max(..mc.map(p => p.canonical_times_t_over_r_over_a))); pairs × probes / $r$ ranges #g(calc.min(..mc.map(p => p.pairs_times_t_over_r))) – #g(calc.max(..mc.map(p => p.pairs_times_t_over_r))). Each arm below is classified by whether its online queries depend on earlier answers (adaptive) or are fixed by the advice and the challenge alone (non-adaptive, DKM's sense: pre-translated queries independent of the challenge).

#text(size: 8pt, table(
  columns: 4, stroke: 0.3pt + luma(190), align: left,
  table.header([arm], [online phase], [class], [what the bound implies]),
  [S3 root Bloom lookup (PR \#1659, `KIC_ROOT_BLOOM=1`)], [Same probe list as the exact root table (fixed state order, keyed by the target's x); the filter only skips exact lookups that would miss; outputs and the 71-probe count are identical in both modes (N23 receipt).], [non-adaptive table scan; Bloom adds ~8 bits per root entry to S, changes no query], [Moves no configuration relative to any bound: S grows by a constant factor, T (probes) is unchanged. At most a constant per probe; it cannot lower the S·T line it sits on.],
  [Multiplier closure (Lead A1)], [Random (c, d) drawn independently of R and of earlier answers; probe T = R − s_c − s_d against the pair table; stop at the first hit.], [non-adaptive], [DKM Thm 1.1: S·T ≳ N below T = √N. Measured canonical·T/(r/A) = #g(calc.min(..mc.map(p => p.canonical_times_t_over_r_over_a))) to #g(calc.max(..mc.map(p => p.canonical_times_t_over_r_over_a))) on all #mc.len() arms: the arm sits on that line and cannot go below it without adaptivity. Matches the earlier "precompute-only" finding (PR \#1648).],
  [Compact-orbit target extraction (n37–n83)], [Fixed state × shift order; S3 roots of the target's x looked up in the K²n-state index; stop at the first lift.], [non-adaptive scan], [At the rank-stage mean, S·T/(r/A) = #g(calc.min(..cr.map(r => r.s_t_over_r_over_a_at_rank_mean))) to #g(calc.max(..cr.map(r => r.s_t_over_r_over_a_at_rank_mean))): on the non-adaptive line. The frozen n71 and n73 targets sit at S·T/(r/A) = #cr.filter(r => r.id.starts-with("compact-n71") or r.id.starts-with("compact-n73")).map(r => g(r.s_t_over_r_over_a)).join(" and "), one draw each, with P(T ≤ frozen) = #cr.filter(r => r.id.starts-with("compact-n71") or r.id.starts-with("compact-n73")).map(r => g(r.p_exponential_at_most_observed)).join(" and ") under an exponential first-hit law at the rank-stage mean.],
  [Pair-root tables / self-collision (Lead B)], [None: relations come from sorting the stored pair states; no target.], [precompute only (no online queries)], [Outside DKM (which bounds online T). Changes the rank-stage precompute only; no measurement exists.],
  [Preprocessing rho (Bernstein–Lange; `KIC_RHO_PRECOMPUTE_WALKS`, `preprocessing_rho.rs`)], [Walk from the target; each step depends on the previous point.], [adaptive], [Not bound by DKM; Corrigan-Gibbs–Kogan S·T² ≈ N is tight. At equal S it needs T ≈ √(N/S), against N/S for any non-adaptive arm. Records are wall-time only (n53, 1,024 targets), so no ops placement here.],
))

== Result 4: multi-target slope

#table(
  columns: 7, stroke: 0.3pt + luma(190), align: (left, left, right, left, right, left, right),
  table.header([source], [curve], [A], [k], [slope], [95% bootstrap], [Kuhn–Struik slope]),
  ..(fit.fits.map(c => ("records", c.curve, str(c.aut), c.ks.map(str).join(","), g(c.slope), "[" + g(c.slope_ci95_bootstrap.at(0)) + ", " + g(c.slope_ci95_bootstrap.at(1)) + "]", g(c.kuhn_struik_law_slope_same_ks))).flatten()
   + meas.fits.map(c => ("fresh", c.curve, str(c.aut), c.ks.map(str).join(","), g(c.slope), "[" + g(c.slope_ci95_bootstrap.at(0)) + ", " + g(c.slope_ci95_bootstrap.at(1)) + "]", g(c.kuhn_struik_law_slope_same_ks))).flatten())
)
Slope = OLS of ln(mean batch total GAE) on ln k; interval = 2,000 bootstrap resamples of batches within each k. Fresh run: all logarithms verified against the planted ones: #meas.fits.map(c => repr(c.all_verified)).join(", ").

== Record checks

#for c in t.record_checks [
  #if c.check == "log2_r" [- #raw(c.where): recorded log2 r = #c.recorded; log2 of the exact r = #g(c.from_registry_r, d: 4).]
  #if c.check == "excluded_precompute" [- #raw(c.where) (n = #c.n): recorded #g(c.recorded_ms) ms excluded; the rung's own components (#raw(c.components_source)): base load #g(c.components_ms.base_load), index build #g(c.components_ms.index_build), rank stage #g(c.components_ms.rank_stage), linear algebra #g(c.components_ms.linear_algebra) ms. The recorded figure equals #c.closest_component_subset.join(" + ") (relative error #g(c.closest_subset_relative_error, d: 5))#if not c.includes_rank_stage [; it leaves out the rank stage].]
]

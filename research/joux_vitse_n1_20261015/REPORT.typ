#set page(paper: "a4", margin: 1.6cm)
#set text(size: 9.2pt)
#set table(stroke: 0.4pt)
#show table: set text(size: 8pt)
#show figure.caption: it => it.body

#let S = json("summary.json")
#let sci(x, d: 2) = {
  if x == none { return "—" }
  let v = float(x)
  if v == 0 { return "0" }
  let e = calc.floor(calc.log(calc.abs(v), base: 10))
  let m = v / calc.pow(10.0, e)
  [#calc.round(m, digits: d)·10#super[#e]]
}
#let r2(x, d: 2) = if x == none { "—" } else { str(calc.round(float(x), digits: d)) }

= Joux–Vitse $(n-1)$-decomposition at $n = 5$: measured against rho

2026-10-15 · `research/joux_vitse_n1_20261015/` · code `src/cryptanalysis/joux_vitse.rs`, `examples/joux_vitse_ladder.rs`, summary by `examples/joux_vitse_summary.rs` · protocol `PROTOCOL.md` (frozen before data) and `PROTOCOL_ADDENDUM.md` · cost model `MODEL.md`. *Class: measurement of a published route (reproduction), with an accounting correction to the scoreboard's $k = 5$ plain constant; no advance.*

== Verdict

No. On $E(FF_(q^5))$ the $(n-1)$-decomposition costs about $12 q^2$ trials of about #sci(S.model_extrapolation.c_trial) $FF_p$ multiplications each (measured, flat in $q$), against rho's $approx 0.89 q^(5\/2)$ group operations of #r2(S.model_extrapolation.c_add, d: 0) multiplications. Measured end to end on complete, verified DLPs it is $4.5 dot 10^7$–$1.8 dot 10^8$ times rho at $q = 11, 31$; on the measured constants the extrapolated crossover is $q^* approx 2^(#r2(S.model_extrapolation.log2_q_star, d: 1))$ (group $approx 2^(#r2(S.model_extrapolation.log2_group_at_q_star, d: 0))$), consistent with Vitse's "faster from a base field larger than $2^60$" (thesis §7.3.3). On the model, IC is still $approx 10^7$ times rho at $q = 2^10$ and $approx 10^6$ at $q = 2^16$ (extrapolation).

== Literature, and two corrections to the premise

The practical $FF_(p^6)$ result of Joux–Vitse (Eurocrypt 2012) is _cover and decomposition_ (a genus-2/3 cover over an intermediate field, Nagao decompositions), not $(n-1)$-decomposition on $E(FF_(p^6))$; the practical $FF_(q^5)$ result of the $n - 1$ variant is an _oracle-assisted static Diffie–Hellman_ attack (Oakley curve, $FF_((2^31)^5)$), not a DLP. The thesis prices one trial at $approx 3 dot 10^9$ multiplications (F4: $10^10$; F4Remake about 3× fewer) and the crossover with rho at $q approx 2^60$ (odd $p$). The repository already had the plain and 2-torsion $k = 5$ test as a stage diagnostic (`jv_quintic.rs`, RESEARCH_K5_TORSION_JOUX_VITSE.md); this study duplicates the plain construction (an incomplete survey, recorded) and adds the trace replay and the first complete $k = 5$ DLPs.

== One trial, across $p$ (#S.rungs.len() rungs, 2 curves × 32 F4 trials each)

#table(columns: 11, align: (right,) * 11,
  [$p$], [$|F|$], [$C_"trial"$ mean], [min], [max], [replay div.], [hits / trials], [trials/rel.], [95 %], [rate·24p], [$C_"rel"$],
  ..S.rungs.map(r => (
    str(r.p), r2(r.base_points_mean, d: 1), sci(r.c_trial_mean), sci(r.c_trial_min), sci(r.c_trial_max),
    str(r.replay_mismatch), str(int(r.rate_hits)) + " / " + str(int(r.rate_trials)),
    r2(r.trials_per_relation, d: 0),
    r2(r.trials_per_relation_wilson95.at(0), d: 0) + "–" + r2(r.trials_per_relation_wilson95.at(1), d: 0),
    r2(r.rate_times_24p), sci(r.c_rel),
  )).flatten()
)

Fits (least squares on logs): $C_"trial" prop p^(#r2(S.fit_c_trial_vs_p.slope))$ (± #r2(S.fit_c_trial_vs_p.slope_se)); trials per relation $prop p^(#r2(S.fit_trials_per_relation_vs_p.slope))$ (± #r2(S.fit_trials_per_relation_vs_p.slope_se)); $C_"rel" prop p^(#r2(S.fit_c_rel_vs_p.slope))$ (± #r2(S.fit_c_rel_vs_p.slope_se)). Checks failed: *#S.checks_failed* (every F4 trial's decompositions equal the meet-in-the-middle set; every relation re-adds to $R$). Replayed trials cost one of three levels ($approx 5.7$, $8.2$, $11.1 dot 10^9$); the suspected cause is `f4_fp` finishing pending pairs as a full run after the recorded steps (not verified per trial). The committed `jv_quintic` plain test re-run on today's engine costs $1.104 dot 10^10$ (it was $3.0 dot 10^10$ on 2026-09-29): the engine's own progress, not this study's.

#figure(image("figures/trial_cost.svg", width: 88%), caption: [Cost of one trial ($FF_p$ multiplications) against $p$; references from Vitse's thesis and PR \#1863's $n = 4$ D#sub[4] oracle (different group).])

#figure(image("figures/trials_per_relation.svg", width: 88%), caption: [Trials per decomposition (Wilson 95 %), the thesis heuristic $24 q$, and the counting floor $p^5 \/ (16 binom(|F|, 4))$ computed from each rung's base.])

== Complete DLPs, paired with rho on the same target

#table(columns: 10, align: (left,) + (right,) * 9,
  [run], [$ell$], [IC trials], [relations], [IC online ($FF_p$ mult.)], [rho ops], [rho $FF_p$ mult.], [IC / rho], [IC / exp. walk], [verified],
  ..S.dlp.map(d => (
    d.file, str(d.ell),
    if d.ic_complete { str(d.ic.trials) } else { "stopped" },
    if d.ic_complete { str(d.ic.relations) } else { str(d.relations_recorded) },
    if d.ic_complete { sci(d.ic.online_muls) } else { "—" },
    str(d.rho_group_ops), sci(d.rho_fp_muls),
    if d.ic_complete { sci(d.ratio_online_ic_over_rho) } else { "—" },
    if d.ic_complete { sci(d.ratio_online_ic_over_expected_rho_walk) } else { "—" },
    if d.ic_complete { if d.ic.verified and d.rho_verified { "both" } else { "NO" } } else { "rho only" },
  )).flatten()
)

Rho's counts include its 20 multipliers' set-up, which dominates at these toy orders (expected walk alone $sqrt(pi ell \/ 2)$); against the walk alone the ratios are larger. $p = 11$, seed 3 drew a curve with $|F| = 1$ and was stopped (`data/dlp/p11_s3.STOPPED.md`). Target-independent set-up (symmetrised $S_5$, $1.87 dot 10^9$; trace recording, $approx 10^10$) is reported separately and excluded from "online".

#figure(image("figures/ic_over_rho.svg", width: 88%), caption: [IC / rho: the complete DLPs (red) and the model (dashed; extrapolation on the measured $C_"trial"$, trials $= kappa q$ with $kappa = #r2(S.model_extrapolation.kappa, d: 1)$, rho with the negation map).])

== $n = 6$, and what was not done

$n = 6$ was not attempted: the symmetrised $S_6$ has 20,349 monomials × 17 powers of $x_R$ (a dense system of $approx 20$ GB against this host's 15 GB, $approx 1.3 dot 10^14$ multiplications), and one trial is six degree-16 equations in five unknowns, Macaulay bound 91 (MODEL.md; derived, a resource limit). No ecbench run: `main` has no $FF_(p^k)$ instance (`Instance::Fp4` is in the unmerged \#1792), recorded in `problems/k5-joux-vitse-trial-cost-and-fpk-ecbench.md`. The curves' ICV1 identities are computed (`curves.json`) but not registered (the registry builder harvests no extension-field JSON).

Environment: Claude Code cloud container, Intel Xeon \@ 2.10 GHz, 4 vCPU, 15 GB, rustc 1.97.0, commit `1e8592da0` for every frozen run; isolation none (wall times are practicality notes only). A container restart killed the first launch; its partial records are kept in `data_killed_20261015T0518Z/` and the restarted runs reproduced their relations exactly.

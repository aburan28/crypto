// Results report. Every number is read from the recorded files under data/;
// nothing here computes a result (rounding is presentation only).
#set page(paper: "a4", margin: (x: 1.5cm, y: 1.4cm), numbering: "1")
#set text(size: 8.8pt)
#set par(justify: true, leading: 0.55em)
#set heading(numbering: "1.")
#show heading: set text(size: 10pt)
#let d = "/data/"
#let q1 = json(d + "random_p2147483647_summary.json")
#let c7 = json(d + "census_p7_summary.json")
#let c11 = json(d + "census_p11_summary.json")
#let ctl = json(d + "control_example1.json")
#let vd = json(d + "verdict.json")
#let ds = json(d + "describe_classes.json")
#let man = json(d + "run_manifest.json")
#let f(x, n: 4) = str(calc.round(x, digits: n))
#let pct(x, n: 2) = str(calc.round(x * 100, digits: n)) + "%"
#let tbl(..args) = text(size: 7.4pt, table(inset: 2.8pt, stroke: 0.4pt + luma(170), ..args))
#let F = q1.fractions

#align(center)[#text(size: 13.5pt, weight: "bold")[Smith's Propositions 2 and 4, natively: per curve and per isogeny class] \
`research/smith_prop24_native_20261014` · 2026-10-14 · protocol `PROTOCOL.md` frozen at commit
#raw(man.at(0).git_head.slice(0, 9)) before any study data · run report `reports/run-20261014T1827Z/`]

= Verdicts (fixed wording of the protocol's decision rules)

#tbl(columns: (20%, 34%, 22%, 24%),
  [*question*], [*measured*], [*rule (frozen)*], [*verdict*],
  [Q1: Smith §8 on random curves, $p = 2^31 - 1$, $N = 10^6$],
  [#F.map(x => f(x.value, n: 4)).join(" / ") (tractable / Prop. 2 / Prop. 4 given Prop. 2 / full_success)],
  [each 98.75% CP interval inside prediction ± 0.01 (#F.map(x => f(x.prediction, n: 4)).join(" / "))],
  [*#q1.overall_verdict*],
  [Q2, $p = 7$ (exhaustive)],
  [$D = #f(c7.direct_full)$, $C = #f(c7.class_full)$, $U = #f(c7.class_tractable)$, $r = #f(c7.r, n: 3)$],
  [$r = (C - D)/(U - D)$: $>= 0.75$ near the ceiling, $<= 0.25$ near the direct fraction],
  [*#c7.verdict*],
  [Q2, $p = 11$ (exhaustive)],
  [$D = #f(c11.direct_full)$, $C = #f(c11.class_full)$, $U = #f(c11.class_tractable)$, $r = #f(c11.r, n: 3)$],
  [same], [*#c11.verdict*],
  [Q2 overall], [both primes], [common verdict, else prime-dependent], [*#vd.q2_overall*],
  [Controls], [Example 1: #if ctl.passed [passed] else [FAILED] · frozen 2026-10-13 counts:
   #if c7.control_frozen_counts_match and c11.control_frozen_counts_match [reproduced] else [NOT reproduced] ·
   class-invariant control: #if c7.control_tractable_closure_equal and c11.control_tractable_closure_equal [closure = direct] else [FAILED] ·
   independence-null median $r$: #f(c7.null.r_median, n: 2) and #f(c11.null.r_median, n: 2) · point count of $X$:
   #(c7.diag.x_agree + c11.diag.x_agree) agree, #(c7.diag.x_disagree + c11.diag.x_disagree) disagree],
  [all pass; instrument invariants all 0], [#if q1.instrument_valid and c7.instrument_valid and c11.instrument_valid [*instrument valid*] else [*invalid*]],
)

*Scope.* This establishes, for $p = 7$ and $p = 11$ only and for hyperelliptic models $y^2 = f(x)$ with $deg f = 8$, the
share of models whose $L$-polynomial class contains a model on which Smith's construction gives an $FF_p$-rational isogeny to
a trigonal-construction Jacobian. It does not establish that any walk reaches that member, what a walk costs, anything
at cryptographic sizes, anything through non-hyperelliptic members of the class, or anything in characteristic 2.
Question 1 reproduces Smith's published rates; it is not new evidence about security.

#figure(image("/figures/smith_prop24_results.svg", width: 96%),
  caption: [(a) Question 1: observed minus predicted, with the recorded 98.75% Clopper–Pearson interval (blue) and
  Smith's Example 2 (orange). (b) Question 2, exhaustive census: direct $D$, class-level $C$ and ceiling $U$; $r$ is
  where $C$ sits between $D$ and $U$. Source `figures/plot.py`, data `data/*_summary.json`.])

= Question 1: Smith §8 reproduced natively

#tbl(columns: 7,
  [*fraction*], [*k / n*], [*value*], [*98.75% CP*], [*cluster bootstrap*], [*prediction (source)*], [*Smith Ex. 2*],
  ..F.map(x => ([#x.name], [#x.k / #x.n], [#f(x.value, n: 5)], [[#f(x.cp_lo, n: 5), #f(x.cp_hi, n: 5)]],
    [#if x.boot_lo == none [—] else [[#f(x.boot_lo, n: 5), #f(x.boot_hi, n: 5)]]],
    [#f(x.prediction, n: 5) (#x.prediction_source)], [#f(x.smith_example2, n: 5)])).flatten())

*Measured.* #q1.n_curves curves carry #q1.subgroups rational tractable subgroups (#f(q1.subgroups_per_curve) per curve; Smith:
1.002244). Every interval contains its prediction and lies inside the ±0.01 margin; Smith's Example 2 value lies in
our interval for #F.filter(x => x.cp_lo <= x.smith_example2 and x.smith_example2 <= x.cp_hi).len() of 4 fractions.
At this $p$ no degenerate case occurred: #q1.diag.maps maps, 0 rank deficiencies, 0 zero discriminants, 0 maps outside
Smith's normal form, 0 degree drops, and the two rational trigonal maps of a subgroup never disagreed on Prop. 4
(#q1.diag.lambda_disagree cases).

= Question 2: per isogeny class at $p = 7$ and $p = 11$

#tbl(columns: 3,
  [*quantity (exact unless marked)*], [*p = 7*], [*p = 11*],
  [models (curve + twist) / $L$-classes], [#c7.curves / #c7.classes], [#c11.curves / #c11.classes],
  [$U$: class has a model with a rational tractable subgroup (= direct)], [#f(c7.class_tractable)], [#f(c11.class_tractable)],
  [$D$: direct `full_success`], [#f(c7.direct_full)], [#f(c11.direct_full)],
  [$C$: class has a `full_success` model], [#f(c7.class_full)], [#f(c11.class_full)],
  [$r = (C-D)/(U-D)$], [#f(c7.r, n: 3)], [#f(c11.r, n: 3)],
  [classes with a tractable model / with a `full_success` model], [#c7.classes_with_tractable / #c7.classes_with_full],
  [#c11.classes_with_tractable / #c11.classes_with_full],
  [independence null: $C$ median [2.5%, 97.5%], 199 permutations], [#f(c7.null.class_fraction_median) [#f(c7.null.class_fraction_lo), #f(c7.null.class_fraction_hi)]],
  [#f(c11.null.class_fraction_median) [#f(c11.null.class_fraction_lo), #f(c11.null.class_fraction_hi)]],
  [measured $P_2$ (Prop. 2 per subgroup) / $P_3$ (Prop. 4 given Prop. 2)], [#f(c7.p2_measured) / #f(c7.p3_measured)],
  [#f(c11.p2_measured) / #f(c11.p3_measured)],
  [Smith's model at this $p$, heuristic $P_2 = (p+1)/(2p)$, $P_3 = 1/2$ (heuristic estimate)],
  [#f(c7.smallp_direct_prediction_heuristic)], [#f(c11.smallp_direct_prediction_heuristic)],
  [Smith's model at this $p$, measured $P_2 P_3$ (heuristic estimate)], [#f(c7.smallp_direct_prediction_measured)],
  [#f(c11.smallp_direct_prediction_measured)],
  [asymptotic: Theorem 2 / $S_8$ ceiling], [#f(c7.asymptotic_direct_prediction) / #f(c7.asymptotic_ceiling)], [same],
)

#let z7 = ds.at(0)
#let z11 = ds.at(1)
*Measured.* At both primes $C$ is far above $D$ and close to $U$: $r = #f(c7.r, n: 3)$ and $#f(c11.r, n: 3)$, both past the
frozen 0.75 threshold, so the verdict is "#vd.q2_overall". $C$ is nevertheless *below* the independence null, which puts
$C$ at $U$ in all 199 permutations at both primes: #z7.tractable_classes_without_full_success of #z7.tractable_classes tractable
classes at $p = 7$ (holding #z7.curves_in_them models) and #z11.tractable_classes_without_full_success of
#z11.tractable_classes at $p = 11$ (#z11.curves_in_them models) contain no `full_success` model, which independent coins would
almost never leave (post hoc, descriptive; Figure 2). Those classes are smaller: median #z7.median_class_size_without_full_success
models against #z7.median_class_size_with_full_success at $p = 7$, #z11.median_class_size_without_full_success against
#z11.median_class_size_with_full_success at $p = 11$. From $p = 7$ to $p = 11$, $r$ rises from #f(c7.r, n: 3) to #f(c11.r, n: 3)
and $D$ moves toward Theorem 2 (#f(c7.direct_full) → #f(c11.direct_full), against 0.1857); Smith's model with the measured
per-subgroup rates gives #f(c7.smallp_direct_prediction_measured) and #f(c11.smallp_direct_prediction_measured) for $D$
(heuristic estimate; consistent with his Hypothesis 4 at these primes).

*Interpretation (labelled as such).* By the frozen rule, the walk extension is worth costing: an isogeny class
that holds a model with a rational tractable subgroup usually also holds a model where both of Smith's square conditions
pass, so a walk inside the hyperelliptic locus could in principle move #pct(c7.class_full) / #pct(c11.class_full) of models to a
`full_success` model, against #pct(c7.direct_full) / #pct(c11.direct_full) directly. The two conditions are not class invariants
(as the tractable subgroup is, by $chi_pi mod 2$), but they are not independent across class members either: a minority of
small classes holds no `full_success` model at all. The reachable ceiling is $C$, not $U$. Nothing here prices the walk
or builds the isogeny between class members.

#figure(image("/figures/class_share.svg", width: 92%),
  caption: [Post hoc, descriptive (protocol amendment 1): share of each tractable class's tractable models that
  have `full_success`, weighted by models. The bin "0" holds the classes with no `full_success` model. Data
  `data/describe_classes.json` from `smith_prop24 describe`.])

= Validation

#tbl(columns: (auto, 1fr),
  [*check*], [*result*],
  [Smith Example 1 ($FF_37$; positive control, same code path)], [#if ctl.passed [passed]: cycle type #ctl.cycle_type.map(str).join("+"),
    #ctl.tractable_subgroups rational tractable subgroup; Smith's $N = x^3 + 16x + 22$, $D = x^2 + 32x + 18$ is one of the
    #ctl.maps.len() computed maps; $f_0, f_1, f_2$, $delta_0, delta_2, delta_4$ and $64 s = delta_1^2$ equal Smith's;
    $\#H(FF_37) = #ctl.h_points$; Prop. 4 rational; the $x -> 1\/x$ model gives the same verdict and its point count of $X$ agrees.],
  [Unit tests (`cargo test --release --lib -- smith_prop24`, 11 tests)], [factorisation against `smith_kernel_coverage`;
    $\#S(H)$ against Lemma 1 and the $S_8$ action; rational lines meeting the four secants counted by brute force over
    $PP^3(FF_7)$ = 2/1/0 by the discriminant's class; $lambda$-roots by a generic determinant; normal-form maps exhausted over
    $FF_7^4$; eq. (7) against the norm; Prop. 4 against the point count of $X$ at $p = 101, 103$; model independence and twist
    flip; Clopper–Pearson reference values.],
  [Point count of Recillas' curve $X$ against Prop. 4 (every census model and twist)],
  [$p = 7$: #c7.diag.x_agree agree, #c7.diag.x_disagree disagree, #c7.diag.x_count_error count errors, #c7.diag.x_not_applicable not applicable;
   $p = 11$: #c11.diag.x_agree / #c11.diag.x_disagree / #c11.diag.x_count_error / #c11.diag.x_not_applicable.],
  [Degenerate cases at small $p$ (handled by the frozen rules)],
  [rank $M < 4$: #c7.diag.m_rank_deficient and #c11.diag.m_rank_deficient subgroups (Prop. 2 fails); zero discriminant:
   #c7.diag.disc_zero and #c11.diag.disc_zero; no degree-3 map although Prop. 2 passes: #c7.diag.construction_failed and
   #c11.diag.construction_failed subgroups (Prop. 4 fails); maps outside Smith's normal form (norm fallback):
   #c7.diag.non_normal_form and #c11.diag.non_normal_form. This is why $P_3 < 1\/2$ at small $p$.],
  [Recomputation], [`smith_prop24 resummarize --dir data` rebuilds every summary byte for byte from the rows.],
)

= Claims by class (controlled-experiment skill §1)

- *Exact computation:* every census count, $D$, $C$, $U$, $r$, the class counts and the Example 1 values.
- *Comparative/statistical:* the four Question 1 equivalences (control: Smith's predictions; 98.75% intervals) and
  the Question 2 verdicts (frozen thresholds on exact values, gated by controls 4 and 5).
- *Post hoc, descriptive (amendment 1):* "$C$ lies below the independence null" (exact $C$ outside the whole range of 199
  stratified permutations) and the within-class spread of Figure 2. Not a registered test.
- *Heuristic estimates:* Theorem 2's 0.1857, the 0.50159 ceiling, and Smith's model at small $p$.
- *Out of scope / not established:* walk cost, explicit isogenies between class members, cryptographic sizes,
  characteristic 2, non-hyperelliptic class members.

= Reproduce

```
cargo build --release --example smith_prop24
B=target/release/examples/smith_prop24; D=research/smith_prop24_native_20261014/data
$B control --out $D
$B random --p 2147483647 --n 1000000 --seed 20261014 --out $D
$B census --p 7 --perms 199 --seed 20261014 --out $D
$B census --p 11 --perms 199 --seed 20261014 --out $D      # about 30 min on 4 vCPU
$B resummarize --dir $D && $B verdict --dir $D && $B describe --dir $D
python research/smith_prop24_native_20261014/figures/plot.py  # layout only
```

*Source.* B. Smith, "Isogenies and the discrete logarithm problem in Jacobians of genus 3 hyperelliptic curves",
J. Cryptology 22 (2009), arXiv:0806.2995, §3, §5, §6, §8. Prior evidence: `research/weak_isogeny_classes_20261013/`.

// Results and run report for research/generic_dlp_correctness_20261017
// (research-visuals skill).  Every number is read from summary.json,
// written by examples/generic_dlp_correctness_summary.rs; all are counts.
#let s = json("summary.json")
#let m = json("raw/manifest.json")
#set page(paper: "a4", margin: 1.7cm, numbering: "1")
#set text(size: 9.5pt)
#set heading(numbering: "1.")
#let f2(x) = str(calc.round(x, digits: 2))

#align(center)[
  #text(15pt, weight: "bold")[Generic-DLP correctness fixes: collab 2-cycle escape, Pohlig–Hellman]

  Study `research/generic_dlp_correctness_20261017` · 2026-10-17 · counts only, nothing timed
]

*What changed.* Findings 3–6 of
`research/generic_dlp_memory_parallelism_20261017/AUDIT_FINDINGS.md`.

+ *`pollard_collab` negation-map 2-cycle escape.*
  - The 2-cycle test compared the unflipped sum with the normalised previous
    point, so it never fired.
  - It now compares after the flip and escapes by doubling the smaller-x
    member of the cycle.
+ *Pohlig–Hellman digits by baby-step giant-step.*
  - Each digit was a linear search over `[0, q)`.
  - Now one `⌈√q⌉`-point table per prime serves all of its digits.
  - Each digit is still the smallest solution below `q`, so the residues
    are unchanged.
+ *Pohlig–Hellman factorisation bounded.* Trial division now stops at the
  smoothness bound.
+ *Pohlig–Hellman answers checked.* `[d]G = Q` is verified, and order 1
  yields `d = 0`.

*Arms.*
- Before: `main` at #raw(m.before.commit).
- After: #raw(m.after.branch).
- Both arms run the same driver, `examples/generic_dlp_correctness_counts.rs`.
- Build: #m.build; #m.rustc.

= Collab walkers

#table(
  columns: 7, stroke: 0.3pt, inset: 3pt,
  [curve], [negation map], [walkers], [DPs before → after], [dead trails before → after], [all DPs verify], [unchanged],
  ..s.collab.map(r => (r.curve, if r.negation_map { "yes" } else { "no" }, str(r.walkers),
    str(r.before.dps) + " → " + str(r.after.dps),
    str(r.before.dead_trails) + " → " + str(r.after.dead_trails),
    if r.all_dps_verified { "yes" } else { "NO" }, if r.unchanged { "yes" } else { "no" })).flatten()
)

- *Without the negation map the code path is not taken*, and every count is
  identical.
- *With it*, walkers that used to die at the step cap in fruitless 2-cycles
  now leave the cycle and reach distinguished points.
- *What still dies* was not classified walker by walker. Longer fruitless
  cycles, which this fix does not handle, are the expected cause, and the
  `20·2^dp_bits` step cap still bounds them.

= Pohlig–Hellman

*The panel.*
- Curve: `y² = x³ + 4` over a 47-bit prime, of order `3·19⁴·151·181·13171`.
- 24 targets × 2 bounds, #s.ph_panel.solves solves in all.
- All correct before and after:
  #if s.ph_panel.all_correct [yes] else [*no*].
- Recovered logarithms identical:
  #if s.ph_panel.recovered_identical [yes] else [*no*].

*Result.* `total_steps` went from #s.ph_panel.total_steps_before to
#s.ph_panel.total_steps_after (#f2(s.ph_panel.ratio)×). Per solve the ratio
ranges #f2(s.ph_panel.per_solve_ratio_min)–#f2(s.ph_panel.per_solve_ratio_max)×.

#figure(image("fig_correctness.svg", width: 100%), caption: [Left: dead trails per curve and mode, before and after. Right: Pohlig–Hellman `total_steps` per solve, before against after (log–log); the dashed line is no change.])

= Pinned walk test re-pin (`tests/ec_vartime_walks_pinned.rs`)

The transcripts the test hashes were captured before and after, in
`transcripts/`. The "before" run used `main`'s sources with the digests
zeroed, and it reproduced the pinned digests exactly.

*Collab.*
- `neg=false` sections byte-identical:
  #if s.pinned.collab_neg_false_identical [yes] else [*no*].
- Dead trails in the `neg=true` sections:
  #s.pinned.collab_neg_true_dead_before → #s.pinned.collab_neg_true_dead_after.

*Pohlig–Hellman and Cheon.*
- Identical apart from step counts:
  #if s.pinned.ph_identical_apart_from_steps [yes] else [*no*].
- The Cheon lines and every answer, factorisation and residue are
  unchanged, wrong-order cases included.
- Over #s.pinned.ph_lines Pohlig–Hellman lines, `total_steps` went from
  #s.pinned.ph_total_steps_before to #s.pinned.ph_total_steps_after.

#table(
  columns: 6, stroke: 0.3pt, inset: 3pt,
  [section], [DP before], [dead before], [DP after], [dead after], [identical],
  ..s.pinned.collab_sections.map(r => (r.section, str(r.dp_before), str(r.dead_before), str(r.dp_after), str(r.dead_after), if r.identical { "yes" } else { "no" })).flatten()
)

= Scope

- *Class.* The collab fix is a correctness fix that changes walk transcripts
  by design. The Pohlig–Hellman change is an algorithmic change to its
  counted cost (`O(e·q)` → `O((e+1)√q)` per prime). Neither touches `S` or
  any IC/rho comparison.
- *What the counts cover.* They cover these curves and targets only. How
  often fruitless 2-cycles arise elsewhere depends on the branch count and
  `dp_bits`.
- *What is not measured.* Pohlig–Hellman wall time; the counts are the
  metric.

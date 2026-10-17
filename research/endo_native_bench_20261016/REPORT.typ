#set page(paper: "a4", margin: (x: 1.8cm, y: 1.8cm), numbering: "1")
#set text(size: 9.5pt)
#set heading(numbering: "1.")
#show table: set text(size: 8pt)

#align(center)[
  #text(size: 15pt, weight: "bold")[Native GLV-chain benchmark: Bandersnatch, Tom-384 and CP6-782 G1]
  #v(2pt)
  Results report · 2026-10-16 · `research/endo_native_bench_20261016/` · harness `85793758e`, protocol and runs `33e2dd10b` · PR \#2006
]

*Summary.* GLV-2 through a cheap isogeny-chain endomorphism, measured natively
with the curve-generic harness `endo-chain-bench` on the three curves the
cross-curve sweep named next (`research/endosweep_curves_20261006`). Counted
field operations ($M_"eq" = M + S$) are the primary result, wall time the
secondary one. All five preregistered conditions of `PROTOCOL.md` are met on
all three curves. Counted ratio best baseline / best GLV-2: *1.525*
(Bandersnatch, short Weierstrass model), *1.521* (Tom-384), *1.492*
(CP6-782 G1), against the frozen model's 1.523 / 1.524 / 1.494. Declared GLV
arm against the fastest baseline by median: 1.504 / 1.504 / 1.570 (minimum
ratios 1.528 / 1.530 / 1.520), A/A median ratios 0.997 / 1.001 / 1.038.
Class *engineering* ($k dot P$ on the same curve, faster); variable-time
harness; no security claim. Every number below is copied from the committed
records and `README.md`; the figures are drawn only from the committed JSON
records by `figures/plot.py`. The run report is
`reports/run-20261016T1759Z/report.pdf`.

= Preregistered conditions

#table(columns: (auto, 1fr, 1fr, 1fr), stroke: 0.4pt + gray,
  [*condition*], [*Bandersnatch*], [*Tom-384*], [*CP6-782 G1*],
  [1. correctness, both runs], [every arm = textbook `BigUint` $k dot P$ on 256/256 pairs; $phi(P) = lambda P$], [yes], [yes],
  [2. $phi(P)$ counts exact; totals within 1 %], [4/4 exact (declared 7M + 2S exception); worst +0.15 %], [44/44; worst −0.24 %], [114/114; worst +0.31 %],
  [3. counted ratio ≥ 99 % of model], [1.525 ≥ 1.508], [1.521 ≥ 1.509], [1.492 ≥ 1.479],
  [4. median > 1.40, timed/counted within ±5 %], [1.504; −1.3 %], [1.504; −1.2 %], [1.570; +4.7 % (inside 2× its 3.8 % A/A floor)],
  [5. lowest $phi(P)$ median = optimised arm at least count], [`0+1w/2`, 226 ns], [`4+1w/7.5.5`, 7 652 ns], [`4+1w/7.5.3`, 28 440 ns],
)

= Counted against the model

#figure(image("figures/fig_counted_vs_model.svg", width: 100%),
  caption: [(a) Counted and model ratio, best baseline / best GLV-2, per curve, with the condition-3 threshold (99 % of model). (b) Deviation of every total arm's counted mean $M_"eq"$ from the frozen model (32, 41 and 48 arms). Source: `<curve>.counts.json`, 256 pairs, seed 20261016.])

The frozen model predicted every total arm from the same recoding and lattice
bases, averaged over 400 other scalars; the counting runs land within 0.31 %
on every arm. The model's counted ratios (1.523 / 1.524 / 1.494) and the
counted ones (1.525 / 1.521 / 1.492) agree to within 0.2 %. The CryptoPro-B
sweep measured 1.408 counted; GLV saves more on the general-$a$ curves, whose
doubling costs 2M + 8S, as the model predicts.

= Timed, with the A/A floor

#figure(image("figures/fig_timed.svg", width: 100%),
  caption: [(a) Declared GLV arm (`optimised/jacobian/w5` on the grid chain) against the reference baseline: counted, median and minimum ratio, hatched band = counted ratio ± the A/A floor $|1 - "A/A median ratio"|$. (b) Per-round ns/op of the reference, the A/A arm and the declared GLV arm, normalised by the reference median (20 rounds each). Source: `<curve>.timing.json`.])

The plotted ratios are against the reference baseline as recorded. On CP6-782
the fastest baseline by median is `jacobian/w5` (2 145 µs), not the reference
`jacobian/w6` (2 175 µs); against it the declared arm's median ratio is 1.570
(condition 4), against the reference 1.592. CP6-782's timing run was the
noisiest: A/A median ratio 1.038, per-round spread 49 % on the reference and
54 % on the A/A arm, 11.1 CPU-s of other activity in the 570 s run (admitted,
`contended: false`, budget 0.5 CPU). Its wall time is read as consistent with
its counts, not as a larger gain.

= Width and table grid

#figure(image("figures/fig_grid.svg", width: 100%),
  caption: [Counted $M_"eq"$ (top) and median µs (bottom) per $k dot P$ against NAF width $w = 3..7$, for the declared chain (optimised evaluator) and the baseline, Jacobian and affine tables. Source: counting and timing JSON.])

Jacobian tables beat affine tables for both arms at every width ≥ 4 on all
three curves; the reason is the Fermat inversion (419, 684 and 1 173 $M_"eq"$).
GLV's optimum is at $w = 5$ everywhere; the baseline's at 5 (Bandersnatch)
and 6 (Tom-384, CP6-782) by count.

= $phi(P)$ stage diagnostics

#figure(image("figures/fig_phi.svg", width: 100%),
  caption: [Every $phi(P)$ arm: counted $M_"eq"$ per call against median µs per call, generic and optimised evaluators (2 + 2, 22 + 22, 57 + 57 chains). Stage diagnostics, never a speed of $k dot P$.])

Every $phi(P)$ arm's counted M and S equals the chain set's model, 162 of 162,
apart from the declared generic `0+1w/2` exception (7M + 2S counted, the
model's 8M + 2S does not correspond to a correct evaluator). Spearman rank
correlation of count and median: 0.949 (4 arms), 0.984 (44), 0.994 (114), as
stated in `README.md`. The chain is 0.2–0.3 %, 2.5 % and 2.2 % of the declared
GLV arm by count.

= Scope and limits

- Bandersnatch is measured on its short Weierstrass model only; its deployed
  twisted-Edwards arithmetic is not measured.
- The baseline is a plain width-$w$ NAF with the same Jacobian arithmetic, not
  a production library; the table conclusion depends on the Fermat inversion.
- Wall-time ratios hold within one curve only; no claim about cycles per field
  operation. Host: Intel Xeon @ 2.10 GHz cloud VM, 4 logical CPUs, timing
  pinned to CPU 3 through `isolated_bench` (L0, not core-isolated hardware).
- Chains outside the frozen chain sets are not measured. Positive control:
  the generic harness reproduces the CryptoPro-B sweep's committed counts on
  all 97 arms with 0 differences (`positive_control/`).

*Records.* `bandersnatch|tom384|cp6_782` × `.counts.json`, `.timing.json`,
`.isolated_bench_record.jsonl`, with their `.stdout.md` tables;
`frozen_inputs/` (chain sets, model, `SHA256SUMS`); `PROTOCOL.md`;
`README.md`. Figures: `figures/plot.py` → `figures/fig_*.svg` / `.png`.

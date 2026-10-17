#set page(paper: "a4", margin: 2cm)
#set text(size: 10pt)
#align(center)[#text(size: 14pt, weight: "bold")[Run report: run-20261017T1623Z-pk-degree-test]\
Observations only. eigenx commit `94fe2d9`; shared 4-vCPU Xeon @ 2.10 GHz VM, rustc 1.97.0, no pinning (L0).]

*Command.* `pk_degree_test L1 100 7 1`, `L3 30 3 0`, `L5 15 3 0`. Per set: one fresh key; the true degree; wrong candidates (uniform, and near misses agreeing on 16, 64 or ⌊f_bits/2⌋ low bits); interleaved timing of one public-key test (theta chain) against one decapsulation. A Lean install and compile overlapped some L1 timing repetitions.

#table(columns: 7, stroke: 0.4pt, inset: 4pt,
  [set], [$E$], [$e$ (true $f$)], [wrong tested / accepted], [test s, median [range]], [decaps s, median [range]], [reps],
  [L1], [1076], [569], [105 / 0], [0.79 [0.53–0.91]], [1.28 [1.11–1.57]], [7],
  [L3], [1594], [832], [35 / 0], [2.72 [2.64–2.78]], [5.54 [5.49–5.75]], [3],
  [L5], [2106], [1087], [18 / 0], [6.84 [6.66–7.03]], [13.34 [11.54–14.26]], [3],
)
Structural inspection at L1: true $f$ glued, 563 Richelot steps, split at the end; one wrong $f$ glued, 566 steps, no early split, no final split. No failures, panics or timeouts; every decapsulation correct.

#figure(image("../../figures/pk_test.svg", width: 95%), caption: [Measured; records in `research/eigen_seta_assessment_20261017/records/run-20261017T1623Z-pk-degree-test/`.])

#set page(paper: "a4", margin: 1.8cm)
#set text(size: 9.5pt)
#set par(justify: true)
#show heading.where(level: 1): set text(size: 14pt)

= ECC2K-130: signed σ steps, $R plus.minus sigma^j (R)$

`ecc2k130/BEST-OF-TWO.md`, branch `research/ecc2k130-best2-walk-20261014`. Runs on Modal, 2026-10-14: emulation `emulation-v{1,2,3}.jsonl` (108 cells, n = 23–59), device reference walk `device-wc/device-wc-v{1,2}.jsonl` (n = 23, 41), RTX PRO 6000 jobs `gpu-run1` (`SIGMA_SIGN=1`, source `f05102952`) and `gpu-run2` (`SIGMA_SIGN=2`). Iteration counts are exact class-collision counts; GPU rates are medians of interleaved 64-launch samples (A/A then A/B), no isolation receipt.

== Result

The campaign σ walk steps $R <- R + sigma^j (R)$. Stepping instead to $R + sigma^j (R)$ or $R - sigma^j (R)$ by a class-invariant sign, the parity of $"popcount"(x and sigma^j x)$ in normal coordinates (`SIGMA_SIGN=2`), changes only the slope's numerator: no extra field product. Measured:

#table(columns: 3, stroke: 0.4pt, inset: 4pt,
  [*quantity*], [*value*], [*source*],
  [rate, RTX PRO 6000, campaign knobs, 385,024 workers], [A/B 0.9965 (A/A drift 0.10%); 14.731 → 14.678 B it/s], [`gpu-run2/results`],
  [iterations, device walk, n = 41], [0.958 (device branches), 0.967 (ECC2K-130 probabilities), ±0.012], [`device-wc-v1/v2`],
  [iterations, emulation, n = 37–59], [0.951–0.959 (ECC2K-130 probabilities)], [`emulation-v3`],
  [correctness], [300/300 replays, 64/64 planted logs, integration + arithmetic pass, 104 regs, 0 spills], [`gpu-run2/results`],
  [GPU time per solve (derived)], [0.961–0.970 of the σ walk], [iterations / rate],
)

Weight bit 4 as the sign (`SIGMA_SIGN=1`) is rate-neutral (A/B 0.99945) but gains $<=$ 1.6% by the model, because the weight is lumpy at n = 131 (minus on 33.9% of steps). No function of the weight can do better. The two walks that *choose* between candidates by a score of the image cut iterations by 13–16%, but they need the second slope: about +38% carry-less work in this kernel (derived, not built).

#figure(image("walk-constant.svg", width: 100%), caption: [Emulated iteration constant $c$ ($plus.minus$ 2 se) against degree, every walk. Left: ECC2K-130 branch probabilities, averaged over mappings. Right: the device's branch choice, one fixed mapping per degree; the σ walk itself moves by several percent between mappings.])

#figure(image("gpu-rate.svg", width: 100%), caption: [RTX PRO 6000 rate samples, A/A panel then A/B panel (red = candidate), for `SIGMA_SIGN=1` (left) and `=2` (right).])

== Caveats

- A different walk: its distinguished points never collide with the current corpus. Switching is the campaign owner's decision, and nothing under the campaign's aws directory is changed.
- The campaign pays for its own n = 131 mapping. Single-mapping rows scatter by several percent, and at n = 53 the signed walk was 6% *worse*. The 3–4% is an average over mappings.
- Short self-returns appear only at n ≤ 37, in the σ walk as well, and none is a formal identity in $ZZ[tau]$.
- `WITNESS=1` is refused with the sign (sixteen counters needed); the campaign runs `WITNESS=0`.

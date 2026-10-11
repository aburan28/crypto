#set page(paper: "a4", margin: (x: 1.8cm, y: 2cm), numbering: "1")
#set text(size: 9.5pt)
#set table(stroke: 0.4pt, inset: 4pt)
#set par(justify: true)
#show heading: set block(above: 1.2em, below: 0.6em)

#align(center)[
  #text(size: 15pt, weight: "bold")[Eigen-Séta v2 in Rust] \
  #text(size: 11pt)[Constant-time decapsulation at NIST level 1, 3 and 5 parameter sizes: results report (two rounds)] \
  #text(size: 8.5pt)[2026-10-11 · runs: run-20261011T1813Z-bench (19c641d12), run-20261011T1827Z-sidechannel (19c641d12), run-20261011T1955Z-verify (d504af28e) · paper: papers/eigen-seta/]
]

= Claims and status

#table(columns: (auto, 2.2fr, 1.8fr, 0.9fr),
  [*\#*], [*Claim*], [*Status*], [*Evidence*],
  [1], [Correct decapsulation at Toy16, Toy32, L1, L3, L5.], [measured: 3,635 decapsulations, no failure], [verify + bench],
  [2], [pk = ct = 1344 / 1992 / 2632 B; decaps 13–15 / 52–54 / 120–124 s at L1 / L3 / L5.], [measured, one host], [bench],
  [3], [The decaps secret phase executes the same instructions for every key and ciphertext, apart from glibc's allocator.], [measured at Toy16 (12 runs) and L1 (4 runs). Evidence, not proof.], [verify, sidechannel],
  [4], [No timing difference between a fixed ct and fresh ones (FO-rejected included), one key.], [measured: \|t\| ≤ 2.36 and ≤ 0.84 in two runs; threshold 4.5], [verify, sidechannel],
  [5], [Robust to adversarial inputs.], [tested: 48 ciphertexts, 14 public keys; no panic or hang], [tests],
  [6], [NIST categories 1 / 3 / 5.], [*conjectured*: binding attack ≈ 2#super[λ] classical, 2#super[λ/2] Grover (shape of AES key search)], [paper §6],
  [7], [The toy `walk` variant is insecure.], [derived from the merged code (key from public degree)], [paper §6.2],
  [8], [The recovered kernel is always cyclic.], [proved (Proposition 3)], [paper §4],
)

Defects found and fixed: an LLVM-introduced conditional in safegcd (found by the callgrind probe); a key-generation hang (Tonelli–Shanks accepted z = 0); degenerate secret scalars ≡ −1 (now rejected at keygen); an uncapped basis search on public curves (now capped). The first version of this report gave quantum security as "about λ/2, as for FESTA" without relating it to NIST's definitions, which are themselves AES key search (claim 6).

= Correctness, sizes and speed

#table(columns: 7,
  [*Set*], [*p bits*], [*pk = ct (B)*], [*keygen (s), n=1*], [*encaps (s)*], [*decaps (s)*], [*correct*],
  [Toy16], [216], [216], [0.38], [0.005–0.008], [0.12–0.19], [5/5],
  [Toy32], [369], [376], [2.3], [0.009–0.012], [0.45–0.59], [3/3],
  [L1], [1341], [1344], [32.3], [0.24–0.26], [13.4–15.0], [3/3],
  [L3], [1991], [1992], [118.6], [1.03], [51.8–53.9], [2/2],
  [L5], [2630], [2632], [271.8], [2.17–2.20], [120.2–124.0], [2/2],
)
The verification run added 3,000/3,000 (Toy16), 600/600 (Toy32) and 20/20 (L1) correct decapsulations with fresh keys. Zero failures in 3,000 trials bounds the failure rate only by about 10#super[−3]. Negligibility rests on the heuristic O(E/p) argument in the paper.

#figure(image("figures/kem_times.svg", width: 78%), caption: [Wall time per operation, all trials, log scale. Measured; run-20261011T1813Z-bench.])

Safegcd inversion is 3.9× to 140× faster than Fermat inversion in the same binary. No KEM-level speedup is claimed, because no baseline binary was kept (AGENTS.md §10).

#figure(image("figures/field_ops.svg", width: 78%), caption: [F_p multiplication and inversions against limb count. Measured; run-20261011T1813Z-bench.])

= Constant-time evidence

#table(columns: 4,
  [*Set*], [*Keys*], [*Runs*], [*Non-allocator instructions*],
  [Toy16], [4], [12], [1,116,209,216 in every run],
  [L1], [2], [4], [120,996,429,610 in every run],
)
Allocator work, including heap growth through sbrk and brk, varies with the heap state left by earlier variable-time phases. It is equal for the two valid ciphertexts of each key.

#figure(image("figures/ct_probe_verify.svg", width: 88%), caption: [Instruction counts of decaps_core per run. Measured; run-20261011T1955Z-verify.])

Limits: instruction counts only, not memory addresses (there are no secret-indexed tables, which is an argument). L3 and L5 were not probed. One compiler and one target. Key generation is variable time by design.

= Key-recovery search space

Key generation's public filters accept 2#super[−8.13], 2#super[−8.64] and 2#super[−9.02] of the odd f candidates at L1, L3 and L5. That leaves about 2#super[268], 2#super[399] and 2#super[525] admissible degrees for an exhaustive search.

#figure(image("figures/f_filter.svg", width: 62%), caption: [log2 of the filter acceptance rate. Measured; run-20261011T1955Z-verify.])

= Graphs checked and left unchanged

docs/index-calculus-scoreboard.html, docs/ic/progress-timeline.json, docs/performance-gains/, docs/curves/ and figures/ were checked. None covers isogeny KEMs. docs/curves/ mentions isogenies only as curve-transfer links in the discrete-log registry, and this work adds no registry curve. This report's graphs are new and are drawn from run records by figures/plot.py and figures/plot_verify.py, which only plot recorded values.

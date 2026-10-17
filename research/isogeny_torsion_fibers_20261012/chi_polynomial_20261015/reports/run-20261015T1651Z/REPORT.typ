#set page(paper: "a4", margin: 19mm)
#set text(size: 10pt)
= χ polynomial paired relation run V1
*Run:* 2026-10-15 16:51 UTC. *Host:* shared Darwin arm64, L0. *Input and source:* `RUN_SPEC_V1.json`.

The enumerative χ polynomial arm solved the same four frozen groups in three rounds each, checked every candidate with exact curve addition, and recovered independently replayed scalars in all 12 runs. Its relation count, rank and trials matched the projected-pair arm in every paired round. Nine raw-pair control slots exhausted. The separate replay checked 99 scalar rows among 108 slots.

#image("GRAPH.svg", width: 100%)

Median whole-pipeline `S_cpu` (χ / projected / negation rho) was 717.54 / 248.45 / 22.04 at `(61,19)`; 293.97 / 142.92 / 12.48 at `(271,61)`; 880.07 / 208.16 / 10.80 at `(787,211)`; and 1441.84 / 214.26 / 10.10 at `(3319,811)`. Paired χ/rho CPU ratios are graphed above; each has three rounds. `S_cpu` is current-thread CPU divided by calibrated source-add CPU and `sqrt(r)`. It is a descriptive CPU equivalent. The reported counted-S field omits the χ pair-equation primitive, so it is not a complete counted method cost.

Raw run SHA-256: `1032c1a29c70d93b18b39dd73c4a305292348f04ed2cbcb39c3df5ca21683ec6`. Independent replay SHA-256: `3aab4f0066ed4bc287dcccda926054f953ca16fa3aa99206a85131dfb6c68d9c`. See `REPORT.md` for exact paths, controls and limits. No isolated wall-time or qualified speed conclusion follows.

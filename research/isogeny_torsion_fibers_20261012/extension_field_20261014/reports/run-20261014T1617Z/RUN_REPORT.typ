#set page(paper: "a4", margin: 16mm)
#set text(size: 9.5pt)
= Run receipt: eight-arm extension-field DLP comparison

*2026-10-14, approximately 16:17 UTC · Darwin 25.6.0 arm64 · shared/L0 · Rust 1.93.1*

The public input is four prime-order subgroups of $E(F_(p^2))$ with orders 19, 61, 211 and 811. The exact fixture SHA-256 is `1a63cc22a22e0fe27d44359325d502038d26d3ac616374fc5558c041939218b7`. The eight-arm V4 spec SHA-256 is `70734d89113dc887a4687560b80f160f1142a5b0cbcd76718190e9b85c01392a`, pinned to release-example source commit `c68835edc4d112e3f8f5da3e606060141ad55143` and source SHA-256 `04ee96dc1af8f40410a550fd36849e54e5d81ea3700c61c81355471118d315e5` before the run.

The V4 raw record SHA-256 is `71a0e8f7d96f630d372299f64a407d398b68dc21b9a8fa04de5c11c8c3b67b6a`. It contains four calibrations and 96 arm rows. The second native arithmetic implementation checked all 96 unique slots and 87 recovered scalars by $[k]G=Q$, retaining nine cap-exhausted raw-pair rows. Its receipt SHA-256 is `59913f5b52e0a54ccb417ee67c57ae8de4aebe91338f470d2aa9a36900c22cb0`. A separate byte-pinned audit confirmed identical relation, query, rank and scalar outcomes among all four projected arms in twelve paired cells.

All exact fields admit the public scaling $psi_u(x,y)=(u^2x,u^3y)$ to the degree-three destination, with $u=27,57,28,264$ and $u^2=-3$. The scaling check passed 64 point-pair homomorphism checks per field plus the frozen generator, target and base images. The V4 run measures this isomorphism arm beside the degree-three arm on the same targets and streams.

The graph reports median cold current-thread CPU divided by a locally calibrated curve-addition time and $sqrt(r)$, with three-round ranges. At order 811 the projected, χ-index, isogeny and scaling medians were 198.17, 277.34, 200.41 and 204.56, against 6.91 for negation rho. The χ/negation-rho paired full-cost ratio was 50.27. The shared host and missing peak RSS prevent an isolated or equal-memory speed claim.

#image("../../GRAPH_V4.svg", width: 100%)

The earlier seven-arm V3 run and its report remain preserved. Exact commands, phase counters, the complete 32-row table, verification scope and open polynomial/pairing obligations are in `RUN_REPORT.md` and `../../REPORT.md`.

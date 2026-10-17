#set page(paper: "a4", margin: 19mm)
#set text(size: 10pt)
= Different-j paired relation solve V1
*Run:* 2026-10-15 17:13 UTC. *Host:* shared Darwin arm64, L0. *Spec:* `RUN_SPEC_V1.json`.

Over $F_(67^2)$ with $s^2=2$, the source $y^2=x^3+46x+37$ has $j=7$ and the explicit degree-three Vélu map reaches $y^2=x^3+53x+25$ with $j=11$. Both have 4575 rational points and order-61 target subgroup. The solver's public input excluded the answer scalar.

All source, isogeny-image and independently built destination factor-base arms attained rank nine and recovered scalar 41 in three of three rounds. Source and image relation outcomes were identical by round: $(22,19)$, $(11,11)$, $(19,17)$ trials and relations. Both negation-rho controls recovered the same scalar. Independent arithmetic checked all 15 arm-round records on both curves.

#image("GRAPH.svg", width: 100%)

The graph shows median complete current-thread CPU in microseconds, three rounds per arm. Source / image / direct destination medians were 151.541 / 178.958 / 216.625; rho / rho A/A were 75.709 / 52.125. Setup, transport, collection, rank solve, recovery and verification are inside each relation timer. The run is L0 and uses a custom native harness without a complete counted-unit conversion.

Source ICV1: `ICV1:fpk-67-2-4083a0ca:-85:4575:7,0:unk:unk:r:653aff0a31ca`. Destination ICV1: `ICV1:fpk-67-2-4083a0ca:-85:4575:11,0:unk:unk:r:3dc527efd80b`. Raw SHA-256: `3e25b3f278977662a33dea7ac932ee49f6de512cec781febd5e43f759c9d0ef2`. See `REPORT.md` for raw and replay links, assumptions and limits.

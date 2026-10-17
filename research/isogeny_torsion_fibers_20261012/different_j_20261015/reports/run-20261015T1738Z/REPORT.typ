#set page(paper: "a4", margin: 19mm)
#set text(size: 10pt)
= Exhaustive rational-fiber audit
*Run:* 2026-10-15 17:38 UTC. *Cell:* $F_(67^2)$, source $j=7$, destination $j=11$, target subgroup order 61. *Spec:* `FIBER_AUDIT_SPEC_V1.json`.

The native audit mapped all 4575 rational source points. All images lay on the destination curve; exactly 1525 distinct image points appeared, each with three source preimages. The map's rational image occupies one third of the 4575-point destination group. All 61 target-subgroup images were distinct and every one of the 3721 ordered subgroup pair-sum checks passed.

#image("GRAPH.svg", width: 100%)

The target subgroup sees an injective transport although the full rational map has three-point fibers. The paired source and isogeny-image relation arms therefore had identical relation outcomes under the frozen projected base. The structural statement is checked on this exact field and map; the full curve-level Lean theorem is still open. `DifferentJ67.lean` proves the stated arithmetic and conditional scalar facts only.

Raw audit SHA-256: `1f193f5ec8bb90f2fb7c742fd6d82cd131410a0fcde688b8f61df5df03d94004`. The recorded audit used 24,142 group additions, 243 doublings, 691,290 field multiplications and 28,387 inversions; these are audit costs, not target-solver costs. See `REPORT.md` for exact source and replay links.

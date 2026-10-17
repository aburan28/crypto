# Direct Semaev root V2 proof map

The [run report](reports/run-20261015T1828Z/REPORT.md) describes an
algorithm and bounded checks. It does not claim a new general summation
polynomial theorem. The readable polynomial identity and exactness scope
are in the [exactness report](reports/run-20261015T1827Z/REPORT.md).

| Mathematical statement | Lean declaration | Exact hypotheses and status |
| --- | --- | --- |
| If an additive image map is injective on its working group, an image relation `image(a)+image(b)=image(t)` is equivalent to `a+b=t`. | `SemaevRootTransfer.imageRelationIff` | Proved for arbitrary types with `Add`, a map-add hypothesis and injectivity. The concrete elliptic curve map is not supplied by this theorem. |
| The four frozen scalar groups of orders 19, 61, 211 and 811 have no nonzero element killed by three. | `kernelVanishes19`, `kernelVanishes61`, `kernelVanishes211`, `kernelVanishes811` | Proved in `Fin r` for exactly the four named orders. |
| If every image collision on one frozen scalar group has difference killed by three, the image is injective. | `injectiveOfKernel19`, `injectiveOfKernel61`, `injectiveOfKernel211`, `injectiveOfKernel811` | Proved conditionally for each order. The kernel-collision premise for the concrete Vélu map remains a separate obligation. |
| The implemented coefficient triple equals the published third summation polynomial and every root corresponds to an unsigned elliptic-curve sum. | Open | The executable coefficients are written in the exactness report and checked by polynomial substitution and exhaustive independent pair enumeration on 22 targets. No Lean polynomial-ring or elliptic-curve model is present in this minimal project; the full universally quantified theorem remains unformalized. |
| The concrete quotient map on all four `F_(p²)` curves satisfies map-add and has the stated kernel on the chosen subgroup. | Open | Native map checks cover the bounded fixture and exactness cases. A faithful Lean curve/field/Vélu definition and kernel proof are still required. |

The project pins Lean `v4.34.1` and has no external dependencies. Build
with `lake build SemaevRootTransfer` in `lean/`. `BUILD_V1.stdout` and
`.stderr` preserve the first failed proof attempt: `decide` on orders
211 and 811 exceeded the recursion limit and reported `sorryAx` for
those failed declarations. V2 replaces those proofs with `omega`;
`BUILD_V2.stdout` records a successful build. Every V2 `#print axioms`
shows either no axioms or only Lean's standard `propext` and
`Quot.sound`; no admitted or custom research axiom appears in the
compiled source. `lean4checker --fresh` was unavailable in this
environment. Neither the V2 build nor the bounded executable checks
close the two concrete curve and polynomial obligations above.

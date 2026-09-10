# Stage 9: degree-21 curve-construction stop

The frozen `K_0` degree-21 public census returned 101 before factor-base construction because `KoblitzCurve::new(0,21)` could not supply the prime-order subgroup required by the current driver. The additive successor changed only the public curve bit to `K_1`; it stopped at the same precondition.

Neither run constructed a factor-base candidate or target, collected a relation, or invoked a solver. The two failed processes consumed 0.002650 and 0.003214 core-seconds, with 0.230542 and 0.009501 seconds wall respectively; both peaked at 1.78 MiB RSS. These are operational failures and support no conclusion about degree-21 factor bases or index calculus.

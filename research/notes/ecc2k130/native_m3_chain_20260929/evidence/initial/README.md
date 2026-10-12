# Initial held attempt: valid but narrow toy witness

The first frozen panel ran from source commit
`16058334e03b139bec49481a8ade541817e707c0` and the accompanying
`FROZEN.json` (copied here), before this archive was evaluated. Its raw
receipt is `PASS`: CryptoMiniSat reported SAT with a complete checked model
for the positive, UNSAT for an exhaustive-oracle negative, and both degree-263
leaf chain representations fit their caps. The independent v1 replay passed.

Post-outcome inspection found that the deterministic first positive target's
solver witness was `T+T+T=T` for the unique x=0 point. It exercised inverse
then copy, while 24 of the 27 frozen toy triples have generic additions at
both edges. Thus this first run does not adequately check the generic-chain
bindings. The original files and hashes are preserved unchanged in this
directory. A separately source-locked revision expands the toy check to all
27 ordered triples and a wrong-target control for each; only that revised
panel may support a semantic PASS decision. The initial result is not a
natural-target yield, n131 solver result, or ECDLP speed claim.

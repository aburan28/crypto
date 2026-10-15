# IC improvement round 2: results (2026-10-15)

These are follow-ups to [`../ic_improve_20261014/README.md`](../ic_improve_20261014/README.md), run under [`RULES.md`](RULES.md):

- Every candidate is compared with the incumbent IC component it replaces.
- Every comparison is paired, uses the same counted unit on both arms, and verifies every logarithm or root set.
- The candidates were tuned on training data, frozen, then tested on held-out data.

| # | Candidate vs incumbent | Verdict | Level | Key numbers | Landed |
|---|---|---|---|---|---|
| [02](02-freerel-production/README.md) | Free relations in the production compact-orbit Koblitz IC, K tuned per arm | **ADVANCES** | constant | <ul><li>Tuned vs tuned, held out: 0.960× [0.946, 0.976] counted, 0.941× wall.</li><li>At n = 61: 0.962× at the same K.</li><li>360/360 runs verified.</li></ul> | **Default on** in `examples/koblitz_orbit_dlp_fast_online.rs`. `KIC_FREE_RELATIONS=0` restores the old pipeline byte-identically (A/A 24/24). |
| [04](04-d4-2lp-fix/README.md) | k = 4 IC: first-dependency stop (`d4-ed`) vs D₄. This fixes round 1's D₄+2LP regression. | **ADVANCES** | constant (exponent MATCHES) | <ul><li>0.979× [0.973, 0.986] over 27 held-out curves, p = 101–1601.</li><li>Never worse on any curve.</li><li>Relations reach the exact first-dependency bound.</li><li>ecbench: 0.990×.</li></ul> | ecbench base `yline-d4-ed`. Existing hashes are unchanged. |
| [05](05-linear-harvest/README.md) | Koblitz Gröbner oracle: harvest the multi-variable linear rows that `decisive_rows` discards, root only, gated on generators per variable < 0.9607 | **ADVANCES** where the gate fires; MATCHES elsewhere | constant | <ul><li>1.075× (19,7) m3; 1.073× (17,9) m2; 1.166× (23,12) m2.</li><li>No cell regresses. 0/96 root-set failures.</li><li>End-to-end reductions −15% / −25%; CPU 1.065× / 0.992×.</li><li>Unit tests 248/248.</li></ul> | **Default on** in `koblitz_groebner.rs`. `KIC_LINEAR_HARVEST=0` restores the old engine (counters 60/60). |
| [03](03-fallrow-gate/README.md) | Koblitz Gröbner oracle: gated fall-row augmentation | ADVANCES where it fires; production no-op | constant | <ul><li>Barrel systems: 1.83× / 4.42× / 11.9× at n = 9 / 11 / 13.</li><li>The gate never fires on pipeline systems, which MATCH the incumbent.</li></ul> | Opt-in `KIC_FALLROW_GATE=1`. |
| [01](01-genus3-linalg/) | Genus-3 2LP IC linear algebra: structured Gaussian elimination, Wiedemann, κ retune | *running* | — | — | — |

## Open follow-ups

- **Koblitz K per size.** Retuning K alone at n = 41 (255 → 100) gives 0.46×, with no new code. The recommended values are K = 125 / 280 / 500 for n = 41 / 53 / 61. Item 02 left the callers' K unchanged; that is an owner decision.
- **Linear-harvest gate.** Generators per variable is close to (N−1)/N on pipeline systems, so the gate skips the largest gain: ungated root harvest gives 1.28× at (29,15) m2. A gate on variable count or tree size needs its own protocol.
- **Linear harvest end to end at (23,12).** The bench gain (1.17×) does not carry through: the production path stops at the first accepted root, so its solves are shorter than the bench's.

## Caveats

- **Host:** shared, 4 vCPU, load average 3–13, isolation level L0. Verdicts use counted operations or thread CPU.
- **Scale:** sizes are small.
- **Conductor:** the control plane was unreachable.

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
| [01](01-genus3-linalg/README.md) | Genus-3 2LP IC: large primes as columns + structured Gaussian elimination + Wiedemann (`lpsgewd`) vs round-1 IC-2LP (DFS substitution + Lanczos on normal equations) | **ADVANCES** | constant; exponent small | <ul><li>Whole solve: 0.653× [0.647, 0.660] at 2^20.</li><li>Exponent q^1.380 vs q^1.398 (Δ −0.018).</li><li>Linear algebra 0.463× at 2^20, but its exponent rises from 1.52 to 1.70 (TRADES).</li><li>The local exponent still rises 0.020 per bit, so it does not settle at 4/3.</li><li>Exploratory κ = 1.41: q^1.357, 0.638× at 2^21.</li><li>474/474 logs verified.</li></ul> | example only |
| [06](06-genus3-kappa-bw/README.md) | Genus-3 2LP IC: frozen κ(q) schedule, block Wiedemann (2,1), and both, vs `lpsgewd` κ = 2 | Block Wiedemann **ADVANCES**; schedule and combination **TRADE** | BW: exponent, constant, primitive; schedule: constant | <ul><li>Block Wiedemann: 0.945× at 2^21, exponent Δ −0.008, linear algebra 0.889× (matrix–vector products 1.50N vs 2.00N).</li><li>Schedule: 0.968× at 2^21 and ≤ 1 at every size, but its fitted exponent rises (+0.028) because the gain is U-shaped. Its local exponent falls to 1.395 at 20–21 (incumbent 1.423), and the trend flattens by 0.011 per bit.</li><li>Combination: 0.929× at 2^21; 0.605× the round-1 IC-2LP.</li><li>2^22 is out of reach: #E overflows u64.</li><li>280/280 logs verified.</li></ul> | example only |

## Open follow-ups

- **Koblitz K per size: done (2026-10-15).** `construct:<n>:<a>` with K omitted now uses the tuned K = 125 / 280 / 500 for n = 41 / 53 / 61 (a = 0). An explicit `:<columns>` still wins, so frozen experiments and CI tune grids are unchanged. Checked by an n = 41 solve: `construct:41:0` and `construct:41:0:125` gave identical summary and target records (timings excluded), 8/8 targets solved. Retuning K alone gives 0.46× at n = 41.
- **Online-only Koblitz rows in `docs/ic/LEADERBOARD.md`: kept** (owner decision, 2026-10-15). They are retired only if there is a strong reason.
- **Linear-harvest gate.** Generators per variable is close to (N−1)/N on pipeline systems, so the gate skips the largest gain: ungated root harvest gives 1.28× at (29,15) m2. A gate on variable count or tree size needs its own protocol.
- **Linear harvest end to end at (23,12).** The bench gain (1.17×) does not carry through: the production path stops at the first accepted root, so its solves are shorter than the bench's.
- **Genus-3 linear algebra.** Elimination removes a shrinking share of columns: it keeps 32% of the incumbent's columns at 2^10 but 91% at 2^20. Weight growth turns into column growth, so linear algebra stays near q^1.54 locally. Item 06 tried the κ retune and block Wiedemann. The relation phase now sets the exponent: its local exponent is about 1.43, still about 0.06 above 4/3.

## Caveats

- **Host:** shared, 4 vCPU, load average 3–13, isolation level L0. Verdicts use counted operations or thread CPU.
- **Scale:** sizes are small.
- **Conductor:** the control plane was unreachable.

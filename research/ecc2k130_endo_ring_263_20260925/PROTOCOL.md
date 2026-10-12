# Protocol: degree-263 kernel direction and orbit cost

Registered 2026-09-25 before running the new certificate. This is a structural and cost-accounting follow-up to merged PRs [#703](https://github.com/aburan28/crypto/pull/703) and [#717](https://github.com/aburan28/crypto/pull/717).

## Frozen inputs

- Exact public curve E: y²+xy=x³+1 over the prescribed F_(2^131), quadratic twist E¹: y²+xy=x³+x²+1, and degree 263.
- Both preflight seeds 20260924 and 20260925, including all eight saved representative generators and the two 263-torsion basis matrices in `research/ecc2k130_direction_review_20260924/twist_torsion_results.json`, SHA-256 `6f1e7b22f3471214d38ec0d3196c88fe9edaf45791e9c05ea67764831fd75368`.
- The independent bit-polynomial field/group arithmetic in `redteam_velu_replay.py`, SHA-256 `c08fbb8de9d54906d783be967a36a6073b2a1d1a3583292badd3adcdd7126790`. It imports no production curve arithmetic.
- The oriented map and field code from #717 for a separate, explicitly implementation-specific operation count. Do not use the oriented map to decide a kernel's direction.

## Hypotheses and exact checks

Derive the degree-two Frobenius polynomial from E(F₂) and compute τ^131=A+Bτ exactly. Check source End(E)=O_(−7), v_263(B)=1, π|E[263]=−I, and the twist's full rational 263-torsion from its exact group order. Evaluate each saved twist generator with independent arithmetic: it has exact order 263, its squared coordinates agree with the saved basis-matrix action, and the line is a τ-eigenline exactly for the two saved horizontal loops. Enumerate all 264 matrix lines per seed; require two fixed lines and two cycles of 131, with the eight representative classifications and codomain b values agreeing with the frozen data. State the codomain order/discriminant for each representative using the standard ordinary-isogeny volcano theorem. Distinguish the group action on the order-r challenge subgroup from a geometric degree-two endomorphism of a fixed leaf.

## Orbit-action accounting

Instrument the #717 map on one fixed public source point for each saved representative. Record field multiplication, squaring, inversion and wall time separately for setup, one full point transport, and one τ=(x²,y²) source action; rerun 11 times after a warmup and report medians and all counts. Keep target validation and codomain membership checks in the map cost. Report the analytic 131-half-kernel loop and the native-leaf alternatives: negation only, generic scalar action, or transport via dual/forward maps with all map/scalar costs charged. These are stage diagnostics, not a matched end-to-end attack comparison.

## Success, failure and stop conditions

Pass only if all exact arithmetic checks agree for both seeds, all eight generators, both matrices and all 264 lines per seed, and operation counters are reproducible. Stop at 180 seconds and preserve a FAIL receipt if any assertion disagrees. Freeze source/input SHA-256, counts, host, Python version, raw medians, and decision in this PR. No factor-base yield, PDP advantage, full ECDLP cost, or n131 speedup may be inferred. The next experiment is a matched native-versus-transported basis/PDP run with cold setup and full rank on the same held-out targets.

## Registered extension before final audit

The certificate also computes the exact minimum norm of a non-scalar element of the index-263 order, with a completing-the-square proof, and verifies the saved generators equal their stated torsion-basis line combinations. On the leaf, verify `phi(tau(P)) = [lambda]phi(P)` at full coordinates and time 11 calls to the existing affine generic scalar routine for `[lambda]`, keeping the scalar routine's field counts separate from map and source Frobenius counts. A deterministic SHA-256 over all exact classification fields and operation counts, excluding wall time and host details, is the CI replay gate. Earlier local runs were developmental; the subsequent frozen receipt is the reported result.

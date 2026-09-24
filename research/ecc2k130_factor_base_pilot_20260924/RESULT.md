# Degree-7 paired factor-base pilot: no native-base advantage on the frozen stream

Status: **toy relation-and-rank diagnostic, negative for this fixed base and stream**. The protocol was committed as `7546c5c` before relation measurements. This experiment uses the degree-7 oriented full-point map from merged [PR #717](https://github.com/aburan28/crypto/pull/717); its exact source SHA-256 is pinned in each raw result. It is not an ECC2K-130 PDP solver or a claim against rho.

The field is `F_(2^21)` with modulus `z^21+z^2+1`. The source has `#E=2,099,948`, trace `-2795`, and an order-421 subgroup coprime to the isogeny degree. Its twist has order `2,094,358`. Deterministic sampling recovered all eight degree-7 twist-kernel lines after 522 abscissae; one returned the source `j` and seven had a different codomain `j`. The predeclared lexicographically first non-self line has kernel abscissae `(142423, 1187869, 1467199)` and codomain `b=0x11584f`. This is a **candidate descending** edge in a ramified degree-7 volcano proxy; it is not a measured degree-263 edge advantage.

Each base has exactly 16 distinct nonzero points after cofactor-4988 projection into the order-421 subgroup. The frozen 512-target stream uses `T_i=[u_i]G+[v_i]Q` with secret `k=356` (present only for audit). The same coefficients are transported to the codomain. An independent paired `F_(q²)` Vélu calculation agreed with the production full-point map on 50 source points, and 16 mapped addition pairs agreed. The standalone verifier independently recomputed all 2,048 variant-target group labels, witnesses, rank stops, base images and recovered scalars with zero discrepancies.

| Base policy | Useful points | Two-summand hits / 512 | First verified rank 17 | Cold field mul to rank | Cold field sqr to rank | Group-add calls to rank |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Original | 16 | 131 | 135 | 89,948 (1.000×) | 92,624 (1.000×) | 4,481 (1.000×) |
| Transported original | 16 | 131 | 135 | 327,008 (3.636×) | 336,153 (3.629×) | 16,573 (3.699×) |
| Descendant native | 16 | 117 | 202 | 360,907 (4.012×) | 371,314 (4.009×) | 18,247 (4.072×) |
| Toy subgroup pullback | 16 | 117 | 202 | 389,771 (4.333×) | 394,745 (4.262×) | 18,668 (4.166×) |

The cold counts include field setup; geometry/kernel search and map construction for all mapped policies; base search and cofactor projection; the exact target-generation prefix to first full rank; one-time subgroup-generator/target transport, then direct codomain target generation where needed; complete pair-table setup; rank work; and scalar verification. The pullback additionally pays for its full 421-point inverse table, an intentionally non-scalable audit device. The field-multiplication count includes operations performed inside inversion; inversion calls are saved separately and **must not** be added again. Group-add counts overlap the field counts and are a second diagnostic unit, not a sum with them. The full 512 attempts, including post-rank cases, remain in the raw JSON. Independent lift and full-stream verification are charged separately as audit phases, not hidden in a candidate cost.

The exact equality of original/transported hit and rank trajectories, and of native/pullback trajectories, is required by the subgroup isomorphism and passed. Merely transporting an existing base supplies no new relations; this pilot spends extra map work for the same rank. The descendant-native base had fewer hits and needed 67 more paired targets to reach full rank. On this one frozen toy stream it fails the predeclared lower-total-cost gate in every reported native unit. That does **not** show that all native bases, other field degrees, or the split degree-263 challenge setting are worse.

The initial raw run is retained as `results_initial.json` under its own historical runner commit `a8e7465`; instrumentation then moved a shared square-root table cost into field setup and added per-target prefix accounting without changing any geometry, base, target, hit, rank, or recovered-scalar record. `results_final.json` is the final charged run. An intermediate local accounting draft double-charged per-target `phi(T)` after forming each source target; the final run corrects this by measuring direct `[u]phi(G)+[v]phi(Q)` generation and treats `phi(T)` as an audit-only equality check. CPU time is in both raw files as a shared-host practicality note; it is not a performance claim. The final raw file is 887,678 bytes with SHA-256 `288b94e14bd4e3b0459160d5c8174f24ee6bbd27928260bbf8d5f9db2d440a23`.

From a fresh checkout of current main:

```sh
python3 research/ecc2k130_factor_base_pilot_20260924/run.py --out /tmp/degree7-new.json
python3 research/ecc2k130_factor_base_pilot_20260924/verify.py /tmp/degree7-new.json
python3 research/ecc2k130_factor_base_pilot_20260924/verify.py research/ecc2k130_factor_base_pilot_20260924/results_final.json
```

The runner refuses to overwrite evidence. Its JSON records source hashes, all 512 cases per variant, phase counters, prefix costs, full-point lift checks, selected kernel data, exact scalar recovery, and limitations. The next useful test would match *multiple* native and original bases by useful cardinality, orbit composition and construction cost on a fresh disjoint target stream, then ask whether the native base reduces cost per new rank. This run supports no base-selection change yet.

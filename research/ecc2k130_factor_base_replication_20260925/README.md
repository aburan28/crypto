# Disjoint-orbit factor-base gate: mixed toy yield, no charged native advantage

[PR #753](https://github.com/aburan28/crypto/pull/753) preregistered this gate at commit `d1e4e56dcd0838f7819a6d38894601cc8a8359a2`, before the new relation measurements. The protocol, two independent base seeds, two disjoint signed-Frobenius target-orbit holdouts, equal 16-point useful subgroup size and rank-17 stop were fixed there. The branch later merged main (including #751) without rewriting that commit. The one-stream [#716](https://github.com/aburan28/crypto/pull/716) result remains a separate, earlier control.

The new degree-7 toy experiment has 16 complete cells: two bases × two 512-accepted-target holdouts × original, transported, descendant-native and pullback policies. All reached rank 17 and recovered the planted scalar 339. Independent code replayed all 8,192 cell outcomes, 2,137 target-candidate draws, four base constructions, pair witnesses, row independence, 16 scalar solutions and the cold cost ledger. Source/transported and native/pullback hit and rank trajectories agree exactly. Source and native base constructions each contain four useful points from each of the same four source signed-τ orbit IDs. Native classification pays for a normalized dual pullback and source-τ canonicalization; it does not pretend that coordinate squaring is a self-map on the leaf.

| Base seed | Holdout | Policy | Hits / 512 | First rank 17 | Cold field mul | / original | Check |
| --- | --- | --- | ---: | ---: | ---: | ---: | --- |
| 2026092511 | A | Original | 153 | 159 | 211,867 | 1.000 | rank, scalar |
| 2026092511 | A | Transported | 153 | 159 | 528,853 | 2.496 | rank, scalar |
| 2026092511 | A | Descendant native | 148 | 97 | 442,132 | 2.087 | rank, scalar |
| 2026092511 | A | Pullback | 148 | 97 | 393,198 | 1.856 | rank, scalar |
| 2026092511 | B | Original | 124 | 148 | 203,749 | 1.000 | rank, scalar |
| 2026092511 | B | Transported | 124 | 148 | 515,565 | 2.530 | rank, scalar |
| 2026092511 | B | Descendant native | 104 | 121 | 490,796 | 2.409 | rank, scalar |
| 2026092511 | B | Pullback | 104 | 121 | 429,542 | 2.108 | rank, scalar |
| 2026092512 | A | Original | 118 | 115 | 167,569 | 1.000 | rank, scalar |
| 2026092512 | A | Transported | 118 | 115 | 462,445 | 2.760 | rank, scalar |
| 2026092512 | A | Descendant native | 132 | 126 | 469,843 | 2.804 | rank, scalar |
| 2026092512 | A | Pullback | 132 | 126 | 406,499 | 2.426 | rank, scalar |
| 2026092512 | B | Original | 198 | 138 | 197,577 | 1.000 | rank, scalar |
| 2026092512 | B | Transported | 198 | 138 | 504,487 | 2.553 | rank, scalar |
| 2026092512 | B | Descendant native | 146 | 289 | 731,005 | 3.700 | rank, scalar |
| 2026092512 | B | Pullback | 146 | 289 | 585,293 | 2.962 | rank, scalar |

The table uses field multiplications, including those inside inversions, as one arithmetic component. Squarings, inversion calls, group-add diagnostics, modular-r operations, CPU and RSS remain separate in `results.json`; group additions overlap field operations and are not added twice. Costs include field and orbit setup, complete kernel search where needed, forward or dual setup, every scanned base candidate/cofactor projection, rejected target draws, target prefix through rank, pair-table construction, missed targets, rank and scalar verification. Native is cheaper than transported in two cells and dearer in two. Its field-multiplication component is 2.087–3.700 times original in all four matched cells, so it fails the preregistered consistent-advantage gate even before an uncalibrated conversion of modular-r work. There is no claim of a calibrated rho or whole-ECDLP speed ratio.

The conditional degree-263 smoke ran because all toy rank and replay gates passed. It uses the two nonconjugate descending lines `[1,0]` and `[1,4]`, the [#750](https://github.com/aburan28/crypto/pull/750) dual, and the actual public ECC2K-130 P,Q from `research/ecc2k130_relations/relations.py` (SHA-256 `0180501ef8fe00f5c54f910822a28cd86dc3c2c1af164348976c1c2e25cc763f`). Their literal coordinates agree with [Certicom's ECC2K-130 parameters](https://www.certicom.com/en/curves-list) and [the independent challenge paper, Appendix A](https://www.ecc-challenge.info/anon.pdf); both satisfy the curve equation and `[r]P=[r]Q=O` in the runner and independent replay. For each line the smoke built 16 useful source and leaf points by cofactor-4 projection, transported/pulled back the controls, charged one native leaf orbit neighbor per selected point via `φ τ [263⁻¹] φ̂`, and checked 32 full-coordinate public `(u,v)` targets. The two lines completed in 24.65 and 24.59 seconds on the recorded, concurrently loaded host, within 300 seconds and 2 GiB. Those wall times are feasibility diagnostics, not comparable performance medians. All four 136-pair m=2 arms had 0/32 hits and rank 0 on both lines; independent code replayed 64 exact targets and 256 arm-target cases. The frozen torsion preflight took 3.267 seconds on its recorded run and is listed separately, since it lacks a compatible field-operation meter. The smoke charges each map and dual construction it performs; no cross-policy cold field-multiplication total is claimed without charging that prior discovery in the same unit.

The m=2 observation is a representation feasibility check. Under uniform subgroup targets, the union-bound expectation is at most `32×136/r = 6.394689268473244e-36` hits per arm. Thus zero hits gives no evidence about an n=131 PDP solver or relation yield. No exact-target scalar was recovered; the n=131 full-DLP cost, `S`, matched rho ratio and speedup remain unset. The degree-7 toy map is ramified, unlike the split degree-263 descent, so its yield trends are not extrapolated.

The next falsifiable gate is an **implicit m≥3 relation producer** with enough supported target mass to make independent useful rank observable before any larger exact leaf run. Preregister original/leaf useful bases, a held-out natural target stream and a counting/memory admission bound; measure S3/PDP misses, transported/native relation parity, tagged-conjugate or actual native leaf-orbit action, cofactor/lift work, rank and same-Q rho. Stop if the admitted base cannot plausibly reach verified rank under the memory cap or if charged native selection remains behind original on both holdouts. The [#751](https://github.com/aburan28/crypto/pull/751) held-out chained affine-screen no-go and [#739](https://github.com/aburan28/crypto/pull/739) n=131 materialized-root bound are hard controls, not bypassed by this smoke.

Reproduce from repository root with fresh outputs:

```sh
base=research/ecc2k130_factor_base_replication_20260925
python3 "$base/run.py" --out /tmp/fb-toy.json
python3 "$base/verify.py" --input /tmp/fb-toy.json --out /tmp/fb-toy-replay.json
python3 "$base/exact_smoke.py" --toy /tmp/fb-toy.json --replay /tmp/fb-toy-replay.json --out /tmp/fb-exact.json
python3 "$base/verify_exact.py" --input /tmp/fb-exact.json --out /tmp/fb-exact-replay.json
python3 "$base/compare_frozen.py" --frozen-toy "$base/results.json" --fresh-toy /tmp/fb-toy.json --frozen-exact "$base/exact_smoke.json" --fresh-exact /tmp/fb-exact.json
```

The committed `results.json`, `replay.json`, `exact_smoke.json`, `exact_replay.json` and `MANIFEST.json` preserve the raw compact outcomes, hashes, phase counters, commands, machine and replay status. Host timings are diagnostic and are excluded from deterministic CI comparison.

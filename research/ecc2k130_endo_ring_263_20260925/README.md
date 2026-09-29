# Degree-263 kernel direction and orbit-action certificate

This independent follow-up to [#703](https://github.com/aburan28/crypto/pull/703) and [#717](https://github.com/aburan28/crypto/pull/717) certifies which of the **actual eight oriented maps tested in #717** descend from the ECC2K-130 crater. The claim is structural. It does not measure PDP yield, relation rank, complete ECDLP cost, or a challenge-scale speedup. The frozen [protocol](PROTOCOL.md), [certificate](certificate.py), and [receipt](receipt.json) make the classification replayable.

For E: `y²+xy=x³+1` over F_(2^131), the degree-two Frobenius `tau` satisfies `tau²+tau+2=0` because E(F₂) has four points. Thus `End(E)=Z[tau]=O_(−7)` is maximal. The exact recurrence gives

```
tau^131 = 8123678690673902378 + 38531015900842053623*tau
38531015900842053623 = 263 * 146505763881528721.
```

The first coefficient is −1 modulo 263 and the second is 0, so F_(2^131)-Frobenius acts as `−I` on all E[263]. A quadratic twist changes this action to `+I`, making its full `(Z/263)²` torsion rational. Its points supply kernel abscissae, but the isogeny still starts on E: a twist point `(x,y)` corresponds over the quadratic extension to `(x,y+s*x)` on E, where `s²+s=1`. Because `s^q=s+1`, the q-Frobenius negates this E point; it preserves every cyclic line. Each line therefore defines an F_q-rational subgroup scheme and quotient even though its nonzero source points are not F_q-rational. The oriented formula in #717 uses only those rational abscissae, and its separate tests check full coordinates and exceptional inputs.

The independent bit-polynomial group implementation from #703 verifies both saved Frobenius matrices on the actual twist basis, all 264 line orbits for each seed, each saved generator's exact 263-order and basis-line identity, and whether squaring its **full point** stays in its own cyclic subgroup. The representative classifications are:

| Preflight seed | Descending representatives | Horizontal controls |
|---|---|---|
| 20260924 | `[1,0]`, `[1,4]` | `[1,39]` (twist eigenvalue 140), `[1,64]` (124) |
| 20260925 | `[1,0]`, `[1,2]` | `[1,30]` (124), `[1,142]` (140) |

In particular, the first saved degree-263 map, line `[1,0]`, is a **descending** map in both seeds; #717 also intentionally covers horizontal controls. All two horizontal lines in each seed have quotient coefficient `b'=1` and return to the source model. The descending coefficients differ from 1 and from one another. The two Frobenius-fixed lines and two 131-cycles per seed match the 264-line classification. The order-theoretic identification uses [Sutherland's ordinary isogeny-volcano theorem](https://arxiv.org/pdf/1208.5370): since 263 splits in O_(−7), exactly two lines preserve the maximal order and the other 262 descend. A descending codomain has `End(E')=Z+263 O_(−7)` (conductor 263, discriminant `−7*263²=−484183`), whereas the source endomorphism discriminant is `−7`. These are different from the much larger discriminant of `Z[tau^131]`.

On a fixed descending leaf, coordinate squaring maps to a distinct Frobenius-conjugate curve, so it is not a degree-two endomorphism of that leaf. In fact, the norm formula `N(a+b*tau)=a²−ab+2b²` with `263|b` gives **121,046** as the least degree of a non-scalar geometric endomorphism of the leaf, attained at `131+263*tau` and `132+263*tau`. This degree bound is not an evaluation-time lower bound. Negation remains cheap. On the public order-`r` subgroup, `tau` acts as the scalar `lambda=196511074115861092422032515080945363956`, satisfying `lambda²+lambda+2=0` and `lambda^131=1 (mod r)`. The certificate independently checks `[lambda]P=tau(P)` and, for every oriented map, `[lambda]phi(P)=phi(tau(P))` at full coordinates. This gives a valid leaf subgroup action, but it needs a charged implementation; it does not restore coordinate squaring on a fixed leaf.

The table reports median stage times over the eight saved maps, each timed 11 times on the recorded macOS/Python host. All field counts include work inside inversion. The point-map setup row starts from a **supplied** kernel generator; torsion discovery and quotient selection are excluded and must be charged in a cold campaign. One source action is measured from 11 batches of 1,024 squarings pairs to limit timer noise.

| Action | Field multiplications | Field squarings | Inversions | Median time |
|---|---:|---:|---:|---:|
| Source `tau(P)` | 0 | 2 | 0 | 2.66 μs |
| Forward oriented `phi(P)` | 1,053 | 400 | 1 | 4.97 ms |
| Leaf generic affine `[lambda]phi(P)` | 25,476 | 25,604 | 193 | 149.08 ms |

The forward map has 131 half-kernel terms. Excluding the 130 multiplies and 131 squarings inside its one Fermat inversion, the measured call uses 923 multiplications and 269 squarings; 917 multiplications and 263 squarings come from the map's main formula and batch inversion, with the remainder in validation. The affine scalar routine is only one available implementation and has not been optimized. These timings do not imply any ECDLP crossover. A source orbit may be labeled before a factor base is transported, but evaluating new native-leaf points or target residuals under the 131-fold action needs scalar, isogeny, or tagged-conjugate work. A horizontal control returns to the source model and retains cheap squaring, so it is not evidence for a new descending representation.

The separate [matched n53 compact-IC/rho comparison in #738](https://github.com/aburan28/crypto/pull/738) and [n131 architecture-specific bound in #739](https://github.com/aburan28/crypto/pull/739) set current attack-cost context; this certificate does not replace either. The next falsifiable test is to freeze the two descending orbit representatives, construct equal-cardinality transported and leaf-native factor bases, and compare held-out PDP success, independent relation rank, and cold full recovery with the orbit-action, target transport, and setup costs charged. A leaf that only looks better when its Frobenius action is treated as free fails that test.

Replay from repository root:

```sh
python3 research/ecc2k130_endo_ring_263_20260925/certificate.py --out /tmp/endo_ring_263_replay.json
```

The frozen receipt has source SHA-256 `6dfe2aaa15b3b819bbee66e4be9afd8d7d737ab71a9cb6fcc13c44cb8f79776f`, exact-input hashes, and deterministic projection SHA-256 `c6129bc8c883244825be7fcfb68a04cfac2c577daf2cd8c39d5063a03c94af0a`. Its local final run passed in 32.74 seconds. Host-specific medians and full per-map counts are retained in the receipt; CI compares the deterministic projection, not timing.

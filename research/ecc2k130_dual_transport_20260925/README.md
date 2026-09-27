# Exact degree-263 dual transport on ECC2K-130

The [protocol](PROTOCOL.md) freezes both saved degree-263 torsion bases from [#703](https://github.com/aburan28/crypto/pull/703). The [research interface](dual_transport.py) constructs a reverse kernel from a complementary twist-torsion point, builds the reverse quotient, and fixes its sign using one non-2-torsion witness. The [certificate](certificate.py) and [final receipt](receipt.json) verify all eight full-point maps tested in [#717](https://github.com/aburan28/crypto/pull/717), including the descending line `[1,0]` in both seeds. [#743](https://github.com/aburan28/crypto/pull/743) independently classifies their directions and explains the leaf endomorphism order.

Let `phi:E→E'` have saved kernel `G=U+sV` on the quadratic twist, where `U,V` are the saved 263-torsion basis. Every tested line has first basis coefficient 1, so `V` lies outside the kernel. Evaluating the same-kernel map on the twist gives `H'=phi_tw(V)`, a rational point of exact order 263 on the twist of `E'`. Its abscissae define `psi:E'→E''`. The normalized reverse quotient has `E''=E` for all eight maps. A separate bit-polynomial full-point implementation reconstructs `H'` through direct rational Vélu sums, rather than trusting the production map to choose its own reverse kernel.

The quotient-kernel argument gives a global identity, not just agreement on sampled points. Under the twist isomorphism, `H'` is the image of the complementary line in `E[263]`. Hence `ker(psi∘phi)=E[263]`, and both `psi∘phi` and `[263]` have degree `263²`. The quotient by this full kernel differs from `[263]` by an automorphism of the ordinary source. Since its nonzero `j` is neither the exceptional `j=0` case nor a larger-automorphism case, that automorphism is `+1` or `−1`. One point whose `[263]` image is not 2-torsion determines the sign. On every saved representative the normalized **raw** reverse map gives `psi∘phi=[−263]`; the formal dual is `[-1]∘psi`, and the certificate checks both `dual∘phi=[263]` and `phi∘dual=[263]` on full public controls. The sign correction is implemented explicitly rather than inferred from x coordinates.

| Saved seed | Descending lines | Horizontal lines | Raw normalized sign | Reverse `b` |
|---|---|---|---:|---:|
| 20260924 | `[1,0]`, `[1,4]` | `[1,39]`, `[1,64]` | −1 on all four | `0x1` |
| 20260925 | `[1,0]`, `[1,2]` | `[1,30]`, `[1,142]` | −1 on all four | `0x1` |

For each of the eight maps, independent direct extension-field Vélu sums agree with the full forward and reverse coordinates of `P`, `Q`, `P+Q`, and `2P`: **64 independent full-coordinate matches**. The dual identity also passes seven signed scalar controls per map, and all four public controls pass the reverse-order identity on `E'`. The exact-target exceptional controls send all 2,096 nonzero forward-twist kernel points and all 2,096 nonzero reverse-twist kernel points to infinity. Infinity passes through all four maps and the formal dual. Sixteen off-curve inputs sharing a kernel abscissa reject before the denominator shortcut, as do eight non-complement and eight sign-ambiguous witness controls. There are no nonzero rational kernel points on the untwisted challenge group; the twist checks exercise the finite exceptions at the actual degree.

The cost measurement uses the existing Python field and oriented-map implementations on the recorded host. It starts with saved `G,V`, so full torsion discovery is **outside** this stage and must be charged in a cold attack. The first row is shared field initialization; construction cost includes all four maps, transport of `V` to the reverse kernel, validation, and the orientation witness. Full composition includes forward transport, raw reverse transport, and sign correction. Medians are across eight representatives; each composition is repeated 11 times after warmup.

| Stage | Field multiplications | Field squarings | Inversions | Median wall |
|---|---:|---:|---:|---:|
| Shared field initialization | — | — | — | 0.268 s |
| Four-map setup and reverse-kernel transport, per representative | 82,638 | 80,750 | 605 | 0.428 s |
| Positive dual composition on one public point | 2,106 | 800 | 2 | 9.09 ms |

Counts include the work inside inversions. Independent replay work is recorded separately and is not a candidate transport cost. The two-map composition count is twice the [#743](https://github.com/aburan28/crypto/pull/743) one-map field count. There is **no** PDP solve, useful relation, factor-base comparison, matrix rank, final scalar recovery, or matched rho comparison here; full-ECDLP cost and speedup remain unset. This settles the exact normalized dual prerequisite, not whether a descended representation helps index calculus.

The initial developmental [FAIL receipt](receipt_initial.json) stopped after 2.66 seconds because the independent checker tested reverse twist points against `b=1` instead of the leaf's `b'`. The construction itself had already produced an order-263 reverse generator and original normalized codomain; correcting the checker made the complete [final receipt](receipt.json) pass in 37.93 seconds. The failure remains visible and is not counted as a mathematical obstruction. Source and input SHA-256 values are in both receipts. The final deterministic projection is `72d6d096b5a0c990b576cfd45ab89973b3396c589ce7165c4e6ff84b891b871f`; CI replays it without comparing host timings.

From repository root, reproduce with a fresh output path:

```sh
python3 research/ecc2k130_dual_transport_20260925/certificate.py --out /tmp/dual_263_replay.json
```

The next falsifiable experiment is the charged, held-out native-versus-transported factor-base/PDP comparison on the two non-conjugate descending representatives. Exact dual transport now makes pullback and round-trip group verification possible; it does not make native-leaf Frobenius free. Charge torsion construction, both maps, leaf orbit action, failed targets, useful rank and final recovery before comparing against matched rho.

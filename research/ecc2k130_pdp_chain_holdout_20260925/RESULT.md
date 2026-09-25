# Chained m=3 PDP affine screen: exact but slower on orbit bases

Status: **toy diagnostic, negative for deployment of this prescreen**. This package follows merged [PR #706](https://github.com/aburan28/crypto/pull/706) and the equal-useful-size degree-7 [PR #716](https://github.com/aburan28/crypto/pull/716). Protocols were committed before their respective target measurements: projected-base bf108d0, raw-base f177f85, masked-base 8ac2c35, and 2,048-target extension 0b44d35. Four compact receipts retain every target outcome and failure, first witness, skip bitset, independent row gain, charged counters, and full-raw SHA-256. Independent replay passed on all four.

The legal placement of the #706 affine contradiction is inside the three-summand chain. For each candidate P3, form U=T−P3; only then apply the two-summand S3(x(P1),x(P2),x(U)) test to the remaining base points. The free chain intermediate is not in the base x-subspace, so applying the original S3 feature to the whole chain would be invalid. The codomain S3 constant is its own b=0x11584f, not 1.

The first projected-subgroup experiment used two disjoint 42-point signed-Frobenius orbit bases and 256 subgroup targets. Each base and its degree-7 image hit **256/256** targets at the first third-point candidate, reached verified rank two at attempt two, and had **zero** certified skips. The tiny order-421 pair table had saturated the target subgroup. The second, raw-point experiment used two disjoint rank-eight x-span orbits and 512 subgroup-only targets; each base and image had **0/512** hits, rank zero, and zero skips. The raw points' cofactor components made subgroup-only targets a poor test of this sparse three-sum set. Both negative runs are preserved.

The third experiment used two morphology-selected, disjoint rank-five raw signed-Frobenius orbits. Each has 42 points and a complete pair table with 903 entries and 757 distinct pair sums. Its complete three-sum set has 7,658 distinct full-group points out of 2,099,948, or 0.3647% exact coverage. A published, independently sampled mask M_i=[421]W_i gives T_i=M_i+[u_i]G+[v_i]Q; [4988]M_i=O preserves the projected-log equation while approximately uniformizing targets over the **full** curve group. The coordinate/sign sampler maps accepted pairs one-to-one to affine points and excludes infinity; that exclusion is the source of its tiny deviation from uniform full-group sampling. This is a synthetic hash-frozen full-group target distribution, not a claim that the ECC2K-130 production target stream already uses torsion masks.

The frozen first 512 masked targets yielded 2/512 hits and verified rank two for the training base at attempt 38, but 0/512 hits and no rank for the held base. The expected count was only 1.87. A separately frozen 2,048-target extension kept the first 512 masks, coordinates, witnesses, skip masks and rank gains identical. Exact coverage predicts 7.47 hits per base at 2,048 approximately uniform targets; the binomial probability of fewer than two is 0.48%. The extension observed 7/2,048 training-base and 4/2,048 held-base group hits. The held base reached independently verified global rank two at attempt 662; its **double-held-out** cell had 2 hits among 649 attempts, giving rank two from that cell alone at global attempt 1871. Independent scalar multiplication verified the recovered toy challenge log k=121 and orbit representative log. Original and transported bases had identical first witnesses, hits, row gains and rank stops on all 2,048 targets.

| 2,048-target policy | Original held base | Transported held base |
| --- | ---: | ---: |
| Confirmed group hits / held target-orbit attempts | 2 / 649 | 2 / 649 |
| Certified residual skips in confirmation | 18,753 | 13,440 |
| Complete pair lookups in confirmation, baseline → screen | 27,193 → 8,440 | 27,193 → 13,753 |
| Field multiplications in confirmation scan, baseline → screen | 601,634 → 1,382,384 (2.30×) | 601,634 → 10,491,134 (17.44×) |
| CPU in confirmation scan, baseline → screen | 0.813 → 3.513 s | 0.819 → 28.194 s |
| Cold field multiplications to first verified global rank, baseline → screen | 1,521,988 → 2,344,853 (1.54×) | 1,295,452 → 11,718,453 (9.05×) |

Cold counts include the 200,000-x morphology search, full pair-table construction, generator/challenge and masked-target generation, and every failed attempt to rank. Transported cold counts also include kernel-line search, map construction and transporting the base and public generator/challenge. The source/transported x-span ranks are 5 versus 20 for the training base and 5 versus 19 for the held base: a point-set isogeny transports relations but does not preserve the low-dimensional abscissa span that makes this screen cheap. Coordinate squaring is not a same-curve automorphism on this codomain, so the transported set is not called a native codomain Frobenius orbit. Curve-addition counts were identical between baseline and screen because both form every residual; the screen saved cheap table lookups while adding field and row-elimination work. CPU is a shared-host diagnostic; field and group counters remain separate native units.

The predeclared confirmation gate **fails** despite an exact screen, no lost hit, and cell-local scalar recovery. Stop promoting the affine S3 screen ahead of a complete pair lookup on these implicit orbit bases. A future test should consider it only where the alternative complete solver is materially more expensive than a lookup, with feature setup, failed targets, rank and rho charged. The degree-7 ramified toy map, order-421 subgroup and 42-point bases cannot establish a degree-263 descendant gain, ECC2K-130 relation coverage, or an n=131/rho crossover.

Replay from a fresh main checkout (Python 3.12):

    base=research/ecc2k130_pdp_chain_holdout_20260925
    python3 "$base/verify.py" "$base/receipt_projected.json"
    python3 "$base/verify.py" "$base/receipt_raw.json"
    python3 "$base/verify.py" "$base/receipt_mask512.json"
    python3 "$base/verify.py" "$base/receipt_power2048.json"
    python3 "$base/verify_prefix.py" "$base/receipt_mask512.json" "$base/receipt_power2048.json"

To regenerate the expensive full raw run without overwriting saved evidence, run python3 "$base/run_mask.py" --attempts 2048 --out /tmp/pdp-power-new.json, then python3 "$base/compact.py" /tmp/pdp-power-new.json /tmp/pdp-power-new-compact.json. The full raw power run took 234.43 s wall on the measurement host; its SHA-256 is dc0bf916453b46d2a17e95eac111d94bbcc133ec18d8db019bc9d69f09942cd6. Each compact receipt embeds code and full-raw hashes; replay independently reconstructs targets, masks, base selection, transport, finite S3 roots, group witnesses, cofactor equations, and rank.

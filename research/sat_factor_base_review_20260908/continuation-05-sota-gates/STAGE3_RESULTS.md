# Stage 3: signed-Frobenius quotient Pollard rho

The rho walk canonicalizes each state over negation and all degree-9 Frobenius images, multiplying both tracked coefficients by the corresponding public factor plus or minus lambda^k. A collision is accepted only when its derived scalar reproduces the public target.

| secret | IC core-s / peak MiB / internal s | automorphism rho core-s / peak MiB / internal s | quotient iterations | walk / reported additions |
|--:|--:|--:|--:|--:|
| 53 | 0.004329 / 2.66 / 0.002152 | 0.003837 / 2.02 / 0.001131 | 2 | 6 / 23 |
| 101 | 0.003912 / 2.58 / 0.002111 | 0.002755 / 2.02 / 0.000826 | 4 | 12 / 29 |

Automorphism rho used fewer total core-seconds in both pairs. Process wall time contains a large unexplained first-launch outlier in each arm, so no wall-time ratio is reported from two samples. Reported additions include the 16 jump-table additions and one initial-state addition; scalar-multiplication internals and Frobenius field operations remain separate costs covered by process timing.

Independent internal review replayed the canonicalization invariant and both collision equations successfully, while returning QUALIFIED on the campaign as a whole. This validates a same-target toy control and the required accounting fields. It does not establish the rho distribution, an exponent, or a cryptographic-size comparison.

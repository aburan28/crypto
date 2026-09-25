# n19 rotated normal-generator support sensitivity

**Decision: the exact sweep admits four fixed beta arms, but none passes the
preregistered follow-up priority gate.** At identical seven-point physical
factors and seven signed cofactor-projected F0 points (three nonzero
negation-paired columns plus O), projected six-sum support ranges from
59,323 to 66,203 of 130,873 subgroup targets. Beta 3 from the merged
[#767 corpus](https://github.com/aburan28/crypto/pull/767) supports 62,389.
The best selected beta gains 3,814 targets (6.11% relative to beta 3), below
the frozen 10% requirement of 68,628 total supported targets. Three arms
improve the exact projected-energy-per-supported-target condition, but the
support conjunct fails for all. The decision is to prioritize other factor-base
and solver questions over more normal-beta tuning on this toy rung.

The [protocol](PROTOCOL.md) was committed at `a6519f6` before selection; the
source, exact four betas and #767 input archive were frozen at `8553b62`
before any candidate six-sum outcome. Its `FROZEN.json` SHA-256 is
`44ac5784e3544c3bc32ed5b63f0c93ff27029c6a12b4598a26160cf492516543`.
Hash-only exact-head CI passed before measurement. The SHA schedule examined
25 counters and selected betas `338435`, `303097`, `464276`, `42605` at
counters `1`, `12`, `13`, `24` respectively. All are normal of rank 19,
with six disjoint rank-two slices, seven physical F0 points, seven signed
projected F0 points and four negation-orbit classes including O. No fallback
was used. The independent bit-serial preflight replay matched all 25 attempts.

Let `N=7^6=117649`, `S` be the exact number of supported projected subgroup
points (including O), and `E=sum_R c_R^2` for complete projected tuple
multiplicities `c_R`. `N^2/E` is the exact Cauchy–Schwarz lower bound on S,
not an observed support size. Pair collisions are `(E-N)/2`; duplicate count
is `N-S`. All five rows have the same necessary counting lower bound of
`q-N=13224` misses, so a change in that bound cannot explain the differences.

| Beta | Distinct full sums | Projected S / q | Exact misses | ΔS vs 3 | Projected duplicate count | E | Colliding tuple pairs | N²/E (approx.) | E/S (approx.) |
|:--|--:|:--|--:|--:|--:|--:|--:|--:|--:|
| 3, frozen #767 reference | 81,175 | 62,389 / 130,873 | 68,484 | — | 55,260 | 327,325 | 104,838 | 42,286.1 | 5.247 |
| 338435 | 81,263 | 66,179 / 130,873 | 64,694 | +3,790 | 51,470 | 297,425 | 89,888 | 46,537.1 | 4.494 |
| 303097 | 82,151 | 66,203 / 130,873 | 64,670 | +3,814 | 51,446 | 298,357 | 90,354 | 46,391.7 | 4.507 |
| 464276 | 82,132 | 59,323 / 130,873 | 71,550 | −3,066 | 58,326 | 353,429 | 117,890 | 39,162.9 | 5.958 |
| 42605 | 80,839 | 64,657 / 130,873 | 66,216 | +2,268 | 52,992 | 309,229 | 95,790 | 44,760.6 | 4.783 |

Every candidate differs by more than the frozen descriptive threshold of
1% of q (1,309 points). This is a structural sensitivity result, not an
algorithmic speed gain. The energy-per-support criterion is checked in exact
integers as `E_candidate*62389 <= 327325*S_candidate`; it passes for
338435, 303097 and 42605 and fails for 464276. None passes the separate
support requirement `S>=68628`, so none passes the conjunctive follow-up gate.

The overlap figures below compare actual supported *point sets* with beta 3,
not independent Bernoulli hit-rate estimates. Candidate-only and beta3-only
sets are disjoint. Common misses equal `q - union`.

| Beta | Intersection with beta 3 | Candidate-only | Beta3-only | Union | Common misses | Fixed #767 labels that flip membership |
|:--|--:|--:|--:|--:|--:|--:|
| 338435 | 31,457 | 34,722 | 30,932 | 97,111 | 33,762 | 5 / 8 |
| 303097 | 31,579 | 34,624 | 30,810 | 97,013 | 33,860 | 5 / 8 |
| 464276 | 28,253 | 31,070 | 34,136 | 93,459 | 37,414 | 4 / 8 |
| 42605 | 30,981 | 33,676 | 31,408 | 96,065 | 34,808 | 4 / 8 |

The eight Q labels were fixed by #767, four planted and four negative for
**beta 3**. For 338435, 303097, 464276 and 42605 respectively, 2, 2, 3
and 2 of those four beta-3 negative labels become supported, while 3, 3, 1
and 2 of its four planted labels become unsupported. The archive contains the
exact `Q`, `R=[4]Q`, all four `Q+T` multiplicities and witnesses for each
candidate, plus its 130,873-entry little-endian target multiplicity array.
Changing the factor base changes these exact memberships; the class names in
the historical #767 target file are not candidate-specific verdicts.

The producer/independent-verifier pair directly covers every labelled tuple
and every subgroup target per beta. Each arm's complete full and projected
histograms, target array, witness identities, overlap and energy were
independently replayed under a separate bit-serial/Fermat field and curve
implementation. All nine children (selection replay plus four producer and
four verifier processes) exited zero, with no cap exceedance or replacement.
Operation types are kept separate: no field-to-point conversion or matched
rho comparator was measured. The candidate producer includes setup, full
census, projection, all-q scan, #767 comparison, witness extraction and raw
serialization; the verifier is additional charged work.

| Beta / child | Point additions | Scalar calls | Field mul | Field square | Field inverse | Wall / CPU s | Peak RSS MiB | Frozen cap |
|:--|--:|--:|--:|--:|--:|:--|--:|:--|
| 338435 producer | 571,591 | 81,317 | 980,544 | 735,472 | 490,235 | 5.027 / 4.920 | 112.1 | 300 s / 512 MiB |
| 338435 verifier | 2,311,501 | — | 17,060,817 | 7,809,689 | 288,464 | 31.288 / 29.984 | 220.0 | 600 s / 512 MiB |
| 303097 producer | 575,339 | 82,205 | 986,256 | 740,988 | 493,091 | 5.247 / 4.989 | 110.0 | 300 s / 512 MiB |
| 303097 verifier | 2,320,615 | — | 17,229,430 | 7,891,718 | 292,254 | 30.379 / 29.732 | 221.3 | 600 s / 512 MiB |
| 464276 producer | 575,263 | 82,186 | 986,144 | 740,876 | 493,035 | 5.471 / 5.166 | 110.4 | 300 s / 512 MiB |
| 464276 verifier | 2,251,655 | — | 16,404,242 | 7,496,344 | 275,436 | 30.328 / 29.272 | 212.3 | 600 s / 512 MiB |
| 42605 producer | 569,727 | 80,893 | 977,664 | 732,760 | 488,795 | 5.364 / 5.165 | 113.2 | 300 s / 512 MiB |
| 42605 verifier | 2,292,041 | — | 16,865,723 | 7,717,003 | 284,712 | 31.429 / 30.368 | 218.2 | 600 s / 512 MiB |

The preflight's own 25-candidate scan used 3,234 field squarings, 330 field
multiplications, 144 inversions, 140 curve additions and 35 scalar calls;
its receipt records 0.007 s wall/CPU and 21.3 MiB peak RSS under the
60 s/128 MiB cap. This wall time is a host-specific diagnostic. Producer
and verifier counts differ because they use separate algorithms and field
implementations; their cost is not combined into an attack-cost number.

The committed [raw archive and receipt](evidence/README.md) include full
source/input/reference hashes, the rejected initial *archive packaging*
attempt and the corrected byte-exact archive. The measured run itself had no
failed child. This experiment measured no solver query, relation rank,
logarithm, end-to-end `S`, rho ratio, n131 support or ECC2K-130 speedup.

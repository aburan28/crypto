#set page(paper: "a4", margin: 1.8cm)
#set text(size: 10pt)
#set par(justify: true)
= NIST seed skip filter: new deterministic features
_2026-10-17 · research/nist_seed_skip_filter_20261017 · preregistration commit 3586810ea_

*Status.* Measured negative result. None of the 21 new preregistered binary features (H7–H11) explains why 11 of the 21 valid NIST/X9.62 seed candidates were not published. This concerns seed provenance, not ECDLP hardness.

*Population.* The 21 candidates accepted by every known criterion in the 2026-10-16 census (`../nist_seed_provenance_20261016`): 10 published, 11 skipped. Features are functions of $(p, n)$, where $t = p+1-n$, $D = t^2 - 4p$ and $n' = 2p+2-n$.

*Hypotheses.* H7: $Z[pi]$ is maximal at every odd $l < 10^6$. H8: the largest prime factor of the twist order (thresholds 160 bits and $"bitlen"(p)-32$; cofactor of at most 32 bits). H9: the same thresholds for $n plus.minus 1$. H10: a shortcut MOV test, $(p/n)$, $(-p/n)$, whether $p$ is a primitive root mod $n$, and whether $(n-1)/"ord"_n(p) <= 2$. H11: $(x/n)$ for $x in {-1,2,3,t,n'}$ and $(y/p)$ for $y in {n,t,D,n'}$.

#figure(image("features.svg", width: 100%), caption: [Recorded values (`out/features.jsonl`, `out/analysis.json`). (a) Largest prime factor of the twist order and of $n plus.minus 1$, in bits relative to $p$; the lines are intervals where factoring did not finish in 600 s. (b) For each binary feature, the share of each group for which it is true. A perfect split would be 1.0 against 0.0.])

*Results (measured).* No binary feature splits the groups perfectly. Strongest leans (partial, not claimed): $(p/n)=1$ in 8/10 published and 5/11 skipped; $"lpf"(n-1) >= 160$ in 2/10 and 6/11. Exact two-sided Mann–Whitney tests over all $binom(21,10)$ labellings: twist $p = 0.51$, $n-1$ $p = 0.20$, $n+1$ $p = 0.52$, count of square factors of $D$ $p = 0.59$. Positive control (planted label) detected. Negative control: 2 of 100,000 random labellings split perfectly on some feature (family-wise rate $2 dot 10^(-5)$). Unfactored after 600 s: P-521 $n+1$, P-384 $n plus.minus 1$; all three are on published curves, which can only block a split. Cross-checks: the largest prime factor of the twist order is 241 bits for P-256 and 118 bits for P-224, both matching known values.

*Caveats.* $(D/p) = 1$ holds by construction ($D equiv t^2$) and was included by mistake. Two H10 conditions are the same condition and are reported once. A bug that inverted the "fully factored" flag was fixed, and the run repeated, before any analysis.

*What remains.* About 80 deterministic features have now been ruled out across both studies. The data fit two hypotheses that exact computation cannot test. One is a randomised check drawn from a generator outside the seed stream (for example, a random-point test $[n]P = O$ without a square-root check, which half of all curves fail). The other is a manual choice. The B-163 seed and four of the base points are still unexplained.

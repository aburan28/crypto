# Exact n19 rotated-base portfolio support

**Decision: pass the preregistered support-only follow-up gate.** The union of
the five fixed equal-size beta arms supports 126,265 of all 130,873 projected
subgroup targets (96.4790%), leaving 4,608 exact misses. The frozen threshold
was `ceil(0.90*q)=117,786`; no three-base subset reaches it, whereas the
fixed-order first four support 121,851 (93.1063%). By the protocol's
fixed-order tie-break, the candidate for a later *charged* four-base study is
`[3,338435,303097,464276]`. The best four-base union in this finite census
is instead `[3,338435,303097,42605]` with 122,355; this is reported as a
descriptive alternative, not substituted into the preregistered decision.

The [protocol](PROTOCOL.md), five-beta order, source and input SHA-256 values
were frozen in commit `a3fd4cb1e749ed939289b2f5b0dbc92c4b889700`
before any three-to-five-base union outcome. Draft PR
[#775](https://github.com/aburan28/crypto/pull/775) and its hash-only CI were
green before the measured run. `FROZEN.json` SHA-256 is
`c5e21c5457ffd66d1d7df752c51ad1a70b3fa466412548289e1b24e0a9012384`;
the two prior raw archives are pinned to
`39f16990213e2c7525c6dae92c131127b12f914bbfd6182264249da72bc9546c`
and `fe84aef6a2cf7f6f4c950245c9c8e870354fb750666f5482b81f9a997d107140`.

| Fixed-order bases admitted | Exact union / q | Coverage | Newly supported by this base | Exact remaining misses |
|:--|--:|--:|--:|--:|
| `3` | 62,389 / 130,873 | 47.6714% | 62,389 | 68,484 |
| `3,338435` | 97,111 / 130,873 | 74.2025% | 34,722 | 33,762 |
| `3,338435,303097` | 114,165 / 130,873 | 87.2334% | 17,054 | 16,708 |
| `3,338435,303097,464276` | 121,851 / 130,873 | 93.1063% | 7,686 | 9,022 |
| all five, then `42605` | 126,265 / 130,873 | 96.4790% | 4,414 | 4,608 |

The previous [single-beta sweep](../rotated_beta_sweep_20260925/RESULT.md)
rejected every *individual* beta under its frozen 10% improvement criterion.
The new result answers a different question: their exact target sets are
complementary. Even beta `464276`, the weakest standalone arm at 59,323
targets, adds 7,686 targets after the first three in the fixed order. Every
one of the 32 five-bit support patterns occurs. The complete pattern counts
and all 31 subset union counts are committed in [outcome.json](evidence/outcome.json).
No Bernoulli independence assumption or sampling interval enters these exact
finite-population counts.

With a *free perfect membership oracle* and the fixed five-base order, a
sequential query would make `258,849/130,873 = 1.977864...` base probes on
average. This is only a number of oracle calls. It does not price constructing
five bases, discovering or certifying support, handling negative PDPs, obtaining
independent logs for the combined columns, relation collection, matrix rank,
or target descent. A solver may pay substantially more on misses; no measured
wall or operation crossover follows from the probe count.

The primary program read complete projected histograms from #767/#769. The
independent bit-serial/Fermat verifier enumerated every `k*[4]H` for
`0<=k<130873`, checked all four #769 `target_counts.u32le` arrays against
their point histograms, mapped the beta-3 histogram through the same group,
then recomputed all 32 patterns and 31 unions. Its comparison of the entire
outcome passed. The cold analysis child took 1.256 s and the independent
verifier 6.274 s on the recorded host; the observed child high-water RSS was
183.9 MiB, below the 120 s / 512 MiB per-child caps. These are archive-analysis
costs, **not** an index-calculus attack cost. The [compact evidence and
receipt](EVIDENCE.md) retain every output and its SHA-256, including
empty stdout/stderr files.

The next experiment should jointly train and solve the four-base logarithm
columns on **new fixed targets** at a rung where complete group-law support
can still certify positives and negatives. Use one target order and compare a
single beta, the frozen four-base choice, and same-Q matched rho; report
base construction, relation yield, independent matrix rank, each negative
PDP attempt, solver setup/search/verification, and independently checked point
logs. A four-base support gate alone is insufficient to start n131 scale-up.
The ongoing n131 F0 exact census and the recursive-S3 exception census must
finish their own semantic gates before such a solver run is admitted. This
experiment measured no rank, solver, ECC2K-130 challenge log, rho ratio,
end-to-end `S`, or n131 transfer.

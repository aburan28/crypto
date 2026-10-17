# N83 K3300 sparse selector plus Bloom prefilter

The sparse canonical selector saved 7.759 ms per target against the
83-rotation scan on 1,000 paired points. A separate root11 Bloom screen saved
3.441 ms per target with the scan selector. This experiment holds the sparse
selector **on in both arms** and toggles only the three-hash 512 MiB Bloom
prefilter. It measures whether Bloom adds an online gain after sparse
canonicalization; it does not add the two separate savings arithmetically.

Both arms use the exact K3300 base of 547,800 usable subgroup points, 3,300
folded columns, base hash
`a4ba33098d2776520b98a326585a2301f8ccd5bc005df7c1a70fbe8496c61d8b`,
903,852,075 exact root keys, fixed five-byte `2^30` root table and
reconstructed 11-byte states. The executable IC source SHA-256 is
`814afa03ae1c386e62e518be3c7a55fb108b09c41b05bef64ca1371787466c39`.
The common public-point file SHA-256 is
`35e34151dcc1ad6048939560adf28fc7c1b1017833fe005a656c3f4e2bc07de8`;
the fixture answer sidecar SHA-256 is
`f8cbc90f7f592aae5014ca4f2693a3327a72e935970f22e8f02868e4aeed7f71`.
These 1,000 points were used in the earlier selector screens, so this is a
secondary paired kernel experiment.

Each arm has a fresh base/index/log process, 14 GiB observed RSS guard and
1,200-second child cap. The timed target interval includes query,
four-summand lookup, relation check, descent and recovery check, while
target-independent setup is recorded separately. The independent verifier
replays the complete base and every scalar and four-point relation and
requires sparse mode in both summaries. The predeclared statistic is the
paired Bloom-on minus Bloom-off online difference over 20 contiguous
50-target blocks, with 10,000 bootstrap resamples and seed `2026101510`.
A negative upper 95% endpoint promotes combined mode for a fresh Tier 1
same-point comparison. Failed, capped and noisy arms remain raw rows.

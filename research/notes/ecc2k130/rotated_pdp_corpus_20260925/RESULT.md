# Exact rotated m5/m6 projected PDP corpus admission

**Decision: PASS for corpus admission only.** The two preregistered arms each
yielded four planted targets and four exact projected-negative targets. The
producer enumerated every labelled tuple; a separate bit-serial field and
group-law implementation directly re-enumerated all 120,774 tuples, rebuilt
both full-point and cofactor-projected histograms, and replayed every target
selection attempt. All four producer/verifier children exited successfully
within their frozen wall and peak-RSS caps. This result freezes point targets
for a later solver comparison; no solver, relation collection, rank computation,
discrete logarithm, or ECC2K-130 attack was measured.

The protocol, source and input hashes were frozen in PR #767 at `39d7327`
before either arm ran, following the merged [rotated-support experiment
#762](https://github.com/aburan28/crypto/pull/762) and [solver-interface
admission #763](https://github.com/aburan28/crypto/pull/763). The frozen
`FROZEN.json` SHA-256 is
`7e0ed24783547ac8a80a83cd62a27976f72798730483121db0c2e55341e45e6b`.
The parent #762 raw archive SHA-256 is
`fad4b4cc48f980349d1ef472d18327368dd0f93f3d1284fc8bbdc20f0d3086d7`.

| Exact arm | q | Points per factor | Labelled tuples | Distinct full sums | Supported projected subgroup points (including O) | Exact misses | Fixed targets |
|:--|--:|:--|--:|--:|--:|--:|:--|
| n13, rotated m5 d2, beta 3 | 2,003 | 5 × 5 | 3,125 | 1,799 | 1,591 / 2,003 (79.43%) | 412 | 4 planted + 4 negative; replay PASS |
| n19, rotated m6 d2, beta 3 | 130,873 | 7 × 6 | 117,649 | 81,175 | 62,389 / 130,873 (47.67%) | 68,484 | 4 planted + 4 negative; replay PASS |

The n19 arm had a preregistered *necessary* lower bound of 13,224 misses:
each factor has at most seven points, so `7^6=117649<q=130873`. Its 68,484
exact misses are a measured group-law census, not a prediction from that bound.
The n13 projected support count agrees point for point with all 2,003 targets
in #762's independently archived parent census; the factor lists also match.
That parent comparison was a post-outcome, read-only crosscheck and did not
change either frozen source or selected targets.

The following is the complete compact target index; the archive retains every
`Q`, projected `R=[4]Q`, planted tuple and torsion shift, all four coset
multiplicities, witnesses, and **every** SHA candidate attempt. For planted
tuples with full sum `S`, the subgroup target is
`Q=[4^(-1) mod q]([4]S)` and the archived shift is `T=S-Q`, with `S=Q+T`.
The verifier checks these identities independently. Each listed negative has
zero full-point multiplicity in all four `Q+T` cosets.

| Arm | Planted SHA counters | Negative SHA counters and subgroup scalars k | SHA attempts |
|:--|:--|:--|--:|
| n13-m5 | 0, 1, 2, 3 | (0, 1764), (4, 274), (6, 588), (13, 806) | 18 |
| n19-m6 | 0, 1, 2, 3 | (0, 83433), (4, 58427), (5, 66830), (7, 17322) | 12 |

All charged operation counts below include setup, factor lifting, the full
histogram and projection, SHA selection, and serialization. The independent
verifier ran separately and is also charged in the receipt. Operation kinds
are reported separately because neither a field-operation conversion nor a
matched rho comparator was measured. Wall and CPU are host-specific diagnostics.

| Arm/child | Point additions | Scalar calls | Field multiplications | Field squares | Field inversions | Wall / CPU s | Peak RSS MiB | Frozen wall/RSS cap |
|:--|--:|--:|--:|--:|--:|:--|--:|:--|
| n13 producer | 10,915 | 1,857 | 18,086 | 15,406 | 9,022 | 0.069 / 0.067 | 23.3 | 300 s / 512 MiB |
| n13 verifier | 46,674 | — | 187,820 | 81,499 | 2,647 | 0.655 / 0.256 | 26.5 | 600 s / 512 MiB |
| n19 producer | 440,796 | 81,221 | 719,158 | 604,779 | 359,537 | 13.982 / 4.513 | 69.3 | 300 s / 512 MiB |
| n19 verifier | 2,142,124 | — | 14,126,741 | 6,417,780 | 224,309 | 35.361 / 27.018 | 203.5 | 600 s / 512 MiB |

The independently enumerated n13 and n19 projected histograms have 1,591 and
62,389 keys, respectively, exactly matching the producer, as do the full
histograms, tuple multiplicity totals, witness identities and all 16 frozen
point labels. No child failure or cap exceedance occurred. The raw archive
SHA-256 is `39f16990213e2c7525c6dae92c131127b12f914bbfd6182264249da72bc9546c`;
its contents, source/input hashes, command, status and replay procedure are
documented in [evidence/README.md](evidence/README.md).

The next gate is a **separately preregistered, byte-for-byte use of these
targets** by candidate and reference solvers, with exact SAT witnesses and
UNSAT verdicts checked by the complete oracle. The present support fractions
do not measure solver difficulty, target-query cost, full ECDLP cost `S`, or a
ratio to rho. In particular, the n19 observed support fraction is a property
of this small curve and these factors, not an n131 or ECC2K-130 forecast.

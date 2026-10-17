# P-192 Type I divisor support and split-base collector

**Result.** The 160 fixed-origin degree-three representatives of the frozen
P-192 Type I signed quotient base have pairwise coprime residual
\(x\)-polynomials: all 12,720 cross-class pairs have gcd one over
\(\mathbb F_p\). The selected residual divisors therefore share no affine
support abscissa. Exactly 20 representatives split into three rational
points over \(\mathbb F_p\). Restricting the collector to those 20 signed
classes leaves 143 geometry-accepted rank-bearing three-sum targets,
compared with 83,874 in the full 160-class base. The frozen 64 independent
accepted known-scalar queries hit none; a steered known-scalar law can hit
64 of 64 known-only targets, but gives rank zero on the original unknown
logarithms.

## Exact divisor geometry

The input classes are the [80 positive](../../type_i_fixed_origin_reduction.json),
[16 negative base](../../type_i_negative_fixed_origin_reduction.json), and
[64 negative guided](../../type_i_guided_negative_fixed_origin_reduction.json)
fixed-origin representatives from the independently checked cubic-section
construction. Each is \([D_i-T]\), with \(D_i\) and the common origin \(T\)
effective rational degree-three divisors on the explicit smooth quartic.
The [Python producer](../../divisor_support_screen.py) normalizes each
residual \(D_i\) as a monic cubic \(h_i(x)\) and a quadratic residual
\(y_i(x)\); it checks squarefreeness, all pairwise polynomial gcds, and
the number of rational roots using \(\gcd(h_i,x^p-x)\). The
[PARI replay](../../verify_divisor_support.gp) independently factors every
cubic and recomputes every gcd from the [frozen GP input](../../type_i_divisor_support_input.gp).
The [receipt verifier](../../verify_divisor_support_receipt.py) checks its
agreement with the original JSON coefficients and the
[GP pass log](../../type_i_divisor_support_gp.log). An initial GP parser
attempt is retained byte-for-byte in its [base64-encoded failed log](../../type_i_divisor_support_gp_attempt1.log.b64);
the corrected replay has the explicit pass marker and no parser errors.

| Signed class group | Classes | Irreducible cubic | One rational root + quadratic | Three rational roots |
| --- | ---: | ---: | ---: | ---: |
| Unknown-log positive | 16 | 8 | 6 | 2 |
| Known-log guided positive | 64 | 26 | 29 | 9 |
| Unknown-log negative | 16 | 6 | 9 | 1 |
| Known-log guided negative | 64 | 19 | 37 | 8 |
| **Total** | **160** | **59** | **81** | **20** |

All 160 cubics are squarefree and monic. The chosen residual supports
contain 141 distinct rational points, 81 degree-two Frobenius orbits, and
59 degree-three orbits, totaling 480 geometric points when orbit degrees
are counted. Pairwise coprimality means none of these residual abscissae
is shared across classes. The common origin \(T\) remains shared by the
representatives; this audit concerns the residual \(D_i\) supports.
The 160 numerator and 160 denominator cubic-section coefficient vectors
are also individually distinct. The [machine receipt](../../type_i_divisor_support.json)
retains each normalized polynomial and splitting result with source and
input hashes.

## Split-only query and relation policy

The [collector](../../split_support_collector.py) selects precisely the
20 classes whose chosen \(h_i\) splits completely. This includes three
signed unknown-log classes and 17 signed known-log classes. It enumerates
all \(\binom{22}{3}=1540\) unordered three-sums with repetition in the
exact P-192 quotient group, records formal unknown-log rows, and applies
the frozen two-character `++` query filter. The [PARI point-law replay](../../verify_split_support_collector.gp)
independently recomputes all 1,540 sums, their known-scalar anchor
equations, their character cells, and the 64 frozen fresh queries. The
[receipt verifier](../../verify_split_support_receipt.py) checks agreement
with the [collector receipt](../../type_i_split_support_collector.json).

| Frozen split-only measurement | Exact value |
| --- | ---: |
| Selected signed classes / unknown-log classes | 20 / 3 |
| Three-sum targets with a nonzero unknown row / known-only targets | 571 / 969 |
| Rank-bearing targets in `++` / known-only targets in `++` | 143 / 282 |
| Repeated target with distinct formal row / bearing-anchor overlap | 0 / 0 |
| Fresh accepted known-scalar queries / rank-bearing hits | 64 / 0 |
| Steered accepted known-scalar queries / raw hits / unknown-row rank | 64 / 64 / 0 |

The 64 steered scalars are derived from known-log triples; PARI checks the
resulting target equals the corresponding scalar multiple of the generator.
The 64 fresh scalars are the disjoint accepted queries from the preceding
[mixed-base run](../run-20261015T1800Z/REPORT.md). For an ideal uniform
accepted nonzero scalar query, the same certified accepted-set lower bound
\(L=1569275433846670190958947353987675545133315652035359680560\)
gives hit probability at most \(143/L<2^{-182}\). This is a support bound
for the stated split-only base and query law. Requiring fully rational
residual support reduces the accepted rank-bearing target set by the exact
factor \(83874/143\), about 586.5, relative to the full signed base.
The cost of constructing and using rational support points is not measured
by this coverage audit, so that reduction is a design tradeoff, not a
runtime comparison.

## Proof map and next collector gate

The pinned [Lean module](../../lean/DivisorSupport.lean) proves that a
polynomial Bézout identity excludes a common root and, when each affine
support lies over a root of its representative polynomial, excludes a
shared affine support point. It also checks the exact integer inequality
\(143\cdot2^{182}<L\). Python and PARI check the concrete P-192 gcd and
splitting instances; the full quartic geometry and 12,720 concrete
Bézout certificates are not encoded as Lean objects. The prior
[`AdaptiveBase` theorem](../../lean/AdaptiveBase.lean) explains why
adjoining private support-point unknowns alone cannot reduce the old-log
dimension. The [run manifest](RUN.json) binds sources, inputs, outputs,
proof logs, and checks by SHA-256.

The next constructive test should reuse a residual support orbit across
independently generated classes, or find a function-divisor relation that
links disjoint supports without introducing a private column for each
class. It then needs ordinary accepted-query yield and projected rank,
followed by target descent and complete one-target cost. The split-only
policy remains a specified collector diagnostic with `candidate_id: null`.

## Reproduction

From the repository root, run the checked producer and independent replays:

```sh
python3 research/p192_extension_cover_20261012/divisor_support_screen.py research/p192_extension_cover_20261012/type_i_fixed_origin_reduction.json research/p192_extension_cover_20261012/type_i_negative_fixed_origin_reduction.json research/p192_extension_cover_20261012/type_i_guided_negative_fixed_origin_reduction.json research/p192_extension_cover_20261012/type_i_divisor_support.json research/p192_extension_cover_20261012/type_i_divisor_support_input.gp
gp -q research/p192_extension_cover_20261012/verify_divisor_support.gp
python3 research/p192_extension_cover_20261012/split_support_collector.py research/p192_extension_cover_20261012/type_i_base_points_legacy_fixture.json research/p192_extension_cover_20261012/type_i_geom_character_collector.json research/p192_extension_cover_20261012/type_i_divisor_support.json research/p192_extension_cover_20261012/type_i_character_density.json research/p192_extension_cover_20261012/type_i_mixed_known_log_support.json research/p192_extension_cover_20261012/type_i_split_support_collector.json research/p192_extension_cover_20261012/type_i_split_support_input.gp
gp -q research/p192_extension_cover_20261012/verify_split_support_collector.gp
```

The two receipt-verifier commands are given by their source usage lines.
Build the pinned proof project with `lake build DivisorSupport`, replay with
`lake env leanchecker --fresh --verbose DivisorSupport`, and audit axioms
with `lake env lean DivisorSupport.lean`.
